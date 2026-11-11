#define NOMINMAX
#include <windows.h>
#include <shellapi.h>
#include <shlobj.h>
#include <shobjidl.h>
#include <objbase.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <sstream>

#include "miniz.h"

namespace fs = std::filesystem;

#ifndef MCPHASE_VERSION
#define MCPHASE_VERSION "6.1"
#endif
static std::wstring utf8_to_wide(const std::string& value);

static std::wstring app_name() {
    return std::wstring(L"McPhase ") + utf8_to_wide(MCPHASE_VERSION);
}




struct Entry {
    std::string destinationUtf8;
    std::uint64_t originalSize = 0;
    std::vector<unsigned char> compressed;
};


static std::wstring utf8_to_wide(const std::string& value) {
    if (value.empty()) return {};

    const int count = MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
        static_cast<int>(value.size()), nullptr, 0);

    if (count <= 0) throw std::runtime_error("Invalid UTF-8 in payload path.");

    std::wstring result(static_cast<std::size_t>(count), L'\0');
    MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
        static_cast<int>(value.size()), result.data(), count);
    return result;
}

static std::uint32_t read_u32(std::istream& in) {
    std::uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        const int c = in.get();
        if (c == EOF) throw std::runtime_error("Unexpected end of payload.");
        value |= static_cast<std::uint32_t>(
                     static_cast<unsigned char>(c)) << (i * 8);
    }
    return value;
}

static std::uint64_t read_u64(std::istream& in) {
    std::uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        const int c = in.get();
        if (c == EOF) throw std::runtime_error("Unexpected end of payload.");
        value |= static_cast<std::uint64_t>(
                     static_cast<unsigned char>(c)) << (i * 8);
    }
    return value;
}

static std::vector<Entry> read_payload(std::istream& in) {


    char magic[4];
    if (!in.read(magic, 4) || std::string(magic, 4) != "MCP1")
        throw std::runtime_error("Invalid payload signature.");

    const std::uint32_t count = read_u32(in);
    if (count > 100000)
        throw std::runtime_error("Payload has an unreasonable file count.");

    std::vector<Entry> entries;
    entries.reserve(count);

    for (std::uint32_t i = 0; i < count; ++i) {
        const std::uint32_t pathSize = read_u32(in);
        const std::uint64_t originalSize = read_u64(in);
        const std::uint64_t compressedSize = read_u64(in);

        if (pathSize == 0 || pathSize > 32768 ||
            compressedSize > 1024ULL * 1024 * 1024 ||
            originalSize > 2ULL * 1024 * 1024 * 1024) {
            throw std::runtime_error("Invalid payload entry size.");
        }

        Entry entry;
        entry.originalSize = originalSize;
        entry.destinationUtf8.resize(pathSize);
        if (!in.read(entry.destinationUtf8.data(), pathSize))
            throw std::runtime_error("Truncated payload path.");

        entry.compressed.resize(static_cast<std::size_t>(compressedSize));
        if (!entry.compressed.empty() &&
            !in.read(reinterpret_cast<char*>(entry.compressed.data()),
                     static_cast<std::streamsize>(entry.compressed.size()))) {
            throw std::runtime_error("Truncated compressed file data.");
        }

        entries.push_back(std::move(entry));
    }

    return entries;
}

static std::vector<Entry> read_embedded_payload() {
    HRSRC resource = FindResourceW(
        nullptr, MAKEINTRESOURCEW(101), MAKEINTRESOURCEW(10));
    if (!resource)
        throw std::runtime_error("Embedded payload resource not found.");

    HGLOBAL loaded = LoadResource(nullptr, resource);
    if (!loaded)
        throw std::runtime_error("Could not load embedded payload.");

    const DWORD size = SizeofResource(nullptr, resource);
    const void* data = LockResource(loaded);
    if (!data || size == 0)
        throw std::runtime_error("Embedded payload is empty.");

    std::string bytes(static_cast<const char*>(data), size);
    std::istringstream stream(
        std::move(bytes), std::ios::in | std::ios::binary);

    return read_payload(stream);
}


static std::vector<unsigned char> decompress_entry(const Entry& entry) {
    if (entry.originalSize > std::numeric_limits<mz_ulong>::max() ||
        entry.compressed.size() > std::numeric_limits<mz_ulong>::max()) {
        throw std::runtime_error("File is too large for this miniz build.");
    }

    std::vector<unsigned char> result(
        static_cast<std::size_t>(entry.originalSize));

    mz_ulong resultSize = static_cast<mz_ulong>(entry.originalSize);
    const int status = mz_uncompress(
        result.data(),
        &resultSize,
        entry.compressed.data(),
        static_cast<mz_ulong>(entry.compressed.size()));

    if (status != MZ_OK || resultSize != entry.originalSize)
        throw std::runtime_error("Could not decompress a payload file.");

    return result;
}

static fs::path safe_destination(const fs::path& root,
                                const std::string& relativeUtf8) {
    std::wstring pathText = utf8_to_wide(relativeUtf8);
    std::replace(pathText.begin(), pathText.end(), L'\\', L'/');

    // DestDir values such as "{app}\\demo" may start with separators
    // after the "{app}" prefix is removed.
    while (!pathText.empty() && pathText.front() == L'/')
        pathText.erase(pathText.begin());

    const fs::path raw(pathText);
    if (raw.is_absolute() || raw.has_root_name() || raw.has_root_directory())
        throw std::runtime_error(
            "Unsafe rooted destination in payload: " + relativeUtf8);

    fs::path clean;
    for (const auto& part : raw) {
        const std::wstring component = part.wstring();

        if (component.empty() || component == L".")
            continue;

        if (component == L".." || component.find(L':') != std::wstring::npos)
            throw std::runtime_error(
                "Unsafe traversal destination in payload: " + relativeUtf8);

        clean /= part;
    }

    if (clean.empty())
        throw std::runtime_error("Empty destination path in payload.");

    return root / clean;
}

static void write_uninstall_manifest(
    const std::vector<Entry>& entries,
    const fs::path& installDir) {

    std::ofstream out(
        installDir / L"uninstall-manifest.txt",
        std::ios::binary | std::ios::trunc);

    if (!out)
        throw std::runtime_error("Cannot create uninstall manifest.");

    for (const auto& entry : entries) {
        // Validate each path before recording it.
        safe_destination(installDir, entry.destinationUtf8);

        out << entry.destinationUtf8 << '\n';
        if (!out)
            throw std::runtime_error("Failed writing uninstall manifest.");
    }
}

static void install_entries(const std::vector<Entry>& entries,
                            const fs::path& installDir) {
    for (const auto& entry : entries) {
        const fs::path target = safe_destination(
            installDir, entry.destinationUtf8);

        fs::create_directories(target.parent_path());
        const auto data = decompress_entry(entry);

        std::ofstream out(target, std::ios::binary);
        if (!out) throw std::runtime_error(
            "Cannot create: " + target.string());

        if (!data.empty()) {
            out.write(reinterpret_cast<const char*>(data.data()),
                      static_cast<std::streamsize>(data.size()));
        }
        if (!out) throw std::runtime_error(
            "Failed writing: " + target.string());
    }
}


static bool create_shortcut(const fs::path& shortcut,
                            const fs::path& target,
                            const fs::path& workingDir,
                            const fs::path& icon) {
    IShellLinkW* link = nullptr;
    HRESULT hr = CoCreateInstance(
        CLSID_ShellLink, nullptr, CLSCTX_INPROC_SERVER,
        IID_IShellLinkW, reinterpret_cast<void**>(&link));

    if (FAILED(hr)) return false;

    link->SetPath(target.c_str());
    link->SetWorkingDirectory(workingDir.c_str());
    if (fs::exists(icon)) link->SetIconLocation(icon.c_str(), 0);

    IPersistFile* persist = nullptr;
    hr = link->QueryInterface(IID_IPersistFile,
                              reinterpret_cast<void**>(&persist));
    if (SUCCEEDED(hr)) {
        hr = persist->Save(shortcut.c_str(), TRUE);
        persist->Release();
    }

    link->Release();
    return SUCCEEDED(hr);
}

static fs::path known_folder(REFKNOWNFOLDERID id) {
    PWSTR raw = nullptr;
    if (FAILED(SHGetKnownFolderPath(id, 0, nullptr, &raw)))
        throw std::runtime_error("Could not locate a Windows folder.");

    fs::path result(raw);
    CoTaskMemFree(raw);
    return result;
}

static fs::path executable_directory() {
    std::wstring buffer(32768, L'\0');
    const DWORD n = GetModuleFileNameW(
        nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));

    if (n == 0 || n >= buffer.size())
        throw std::runtime_error("Could not locate the installer executable.");

    buffer.resize(n);
    return fs::path(buffer).parent_path();
}

static std::string lower_ascii(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });
    return s;
}

static void show_version(const std::vector<Entry>& entries) {
    for (const auto& entry : entries) {
        std::string path = lower_ascii(entry.destinationUtf8);
        std::replace(path.begin(), path.end(), '\\', '/');

        if (path == "version" ||
    (path.size() >= 8 &&
     path.compare(path.size() - 8, 8, "/version") == 0)) {

            const auto bytes = decompress_entry(entry);
            std::cout << "Installer information:\n";
            if (!bytes.empty())
                std::cout.write(reinterpret_cast<const char*>(bytes.data()),
                                static_cast<std::streamsize>(bytes.size()));
            std::cout << "\n\n";
            return;
        }
    }
}

int main() {
    const HRESULT comResult = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);

    try {
        const auto entries = read_embedded_payload();
        show_version(entries);

        std::wstring version = utf8_to_wide(MCPHASE_VERSION);
std::replace(version.begin(), version.end(), L'.', L'_');
std::wstring installText = L"C:\\mcphas" + version;
std::wcout << L"Install directory [" << installText << L"]? ";
        std::wstring input;
        std::getline(std::wcin, input);
        if (!input.empty()) installText = input;

        const fs::path installDir(installText);
        std::wcout << L"Create desktop shortcuts? [y] ";
        std::getline(std::wcin, input);
        const bool desktopIcons =
            input.empty() || (input[0] == L'y' || input[0] == L'Y');

        std::wcout << L"Installing to " << installDir << L"...\n";
        std::wcout << L"Extracting files...\n";
        install_entries(entries, installDir);
         std::wcout << L"... complete.\n";
        std::wcout << L"Creating Uninstaller ...\n";

write_uninstall_manifest(entries, installDir);

HRSRC uninstallerResource = FindResourceW(
    nullptr, MAKEINTRESOURCEW(102), MAKEINTRESOURCEW(10));

if (!uninstallerResource)
    throw std::runtime_error("Embedded uninstaller resource not found.");

HGLOBAL uninstallerLoaded = LoadResource(nullptr, uninstallerResource);
const DWORD uninstallerSize = SizeofResource(nullptr, uninstallerResource);
const void* uninstallerData =
    uninstallerLoaded ? LockResource(uninstallerLoaded) : nullptr;

if (!uninstallerData || uninstallerSize == 0)
    throw std::runtime_error("Embedded uninstaller is empty.");

std::ofstream uninstallerOut(
    installDir / L"uninstall.exe",
    std::ios::binary | std::ios::trunc);

if (!uninstallerOut)
    throw std::runtime_error("Cannot create uninstall.exe.");

uninstallerOut.write(
    static_cast<const char*>(uninstallerData),
    static_cast<std::streamsize>(uninstallerSize));

if (!uninstallerOut)
    throw std::runtime_error("Failed writing uninstall.exe.");

          std::wcout << L"... complete.\n";
        std::wcout << L"Creating Shortcuts ...\n";


        const std::wstring appName = app_name();

const fs::path programs =
    known_folder(FOLDERID_Programs) / appName;
fs::create_directories(programs);

create_shortcut(
    programs / (appName + L".lnk"),
    installDir / L"mcdos.bat",
    installDir,
    installDir / L"mcphas.ico");

create_shortcut(
    programs / (appName + L" - Demo.lnk"),
    installDir / L"demo" / L"demo.bat",
    installDir / L"demo",
    installDir / L"demo" / L"demo.ico");

create_shortcut(
    programs / (appName + L" - Uninstall.lnk"),
    installDir / L"uninstall.exe",
    installDir,
    installDir / L"mcphas.ico");

  std::wcout << L"... complete.\n";

        if (desktopIcons) {
        std::wcout << L"Creating Desktop Icons ...\n";
            const fs::path desktop = known_folder(FOLDERID_Desktop);

            create_shortcut(
                desktop / (appName + L".lnk"),
                installDir / L"mcdos.bat",
                installDir,
                installDir / L"mcphas.ico");

            create_shortcut(
                desktop / (appName + L" - Demo.lnk"),
                installDir / L"demo" / L"demo.bat",
                installDir / L"demo",
                installDir / L"demo" / L"demo.ico");
     std::wcout << L"... complete.\n";
   }

        std::wcout << L"Setting up startup by running makedos.bat ... ";
          const std::wstring params = installDir.wstring();
            ShellExecuteW(nullptr, L"open",
                          (installDir / L"makedos.bat").c_str(),
                          params.c_str(), installDir.c_str(), SW_SHOWNORMAL);
           std::wcout << L"... complete.\n";
  

        std::wcout << L"Install McPhase screensaver? [y] ";
        std::getline(std::wcin, input);
        if (input.empty() || (input[0] == L'y' || input[0] == L'Y')) {
          std::wcout << L"Installing McPhase Screensaver ...\n";
            const fs::path screensaver =
                installDir / L"demo" / L"McPhase_screensaver.scr";
            const std::wstring params =
                L"desk.cpl,InstallScreenSaver \"" +
                screensaver.wstring() + L"\"";
            ShellExecuteW(nullptr, L"open", L"rundll32.exe",
                          params.c_str(), nullptr, SW_SHOWNORMAL);
           std::wcout << L"... complete.\n";
  
        }
    } catch (const std::exception& e) {
        std::cerr << "Installation failed: " << e.what() << '\n';
        if (comResult == S_OK || comResult == S_FALSE) CoUninitialize();
        return 1;
    }
  
    if (comResult == S_OK || comResult == S_FALSE) CoUninitialize();
    std::wcout << L"McPhase Installation complete.\n";

    return 0;
}
