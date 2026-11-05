#define NOMINMAX
#include <windows.h>
#include <shlobj.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>

namespace fs = std::filesystem;

#ifndef MCPHASE_VERSION
#define MCPHASE_VERSION "6.1"
#endif

static fs::path executable_path() {
    std::wstring buffer(32768, L'\0');
    const DWORD n = GetModuleFileNameW(
        nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));

    if (n == 0 || n >= buffer.size())
        throw std::runtime_error("Could not locate the uninstaller.");

    buffer.resize(n);
    return fs::path(buffer);
}

static fs::path known_folder(REFKNOWNFOLDERID id) {
    PWSTR raw = nullptr;
    if (FAILED(SHGetKnownFolderPath(id, 0, nullptr, &raw)))
        throw std::runtime_error("Could not locate a Windows folder.");

    fs::path result(raw);
    CoTaskMemFree(raw);
    return result;
}

static std::wstring app_name() {
    const std::string utf8 = MCPHASE_VERSION;
    const int count = MultiByteToWideChar(
        CP_UTF8, 0, utf8.c_str(), -1, nullptr, 0);

    if (count <= 1)
        throw std::runtime_error("Invalid version string.");

    std::wstring version(static_cast<size_t>(count), L'\0');
    MultiByteToWideChar(
        CP_UTF8, 0, utf8.c_str(), -1, version.data(), count);
    version.resize(static_cast<size_t>(count - 1));

    return L"McPhase " + version;
}


static fs::path safe_destination(
    const fs::path& root, const std::string& relativeUtf8) {

    if (relativeUtf8.empty())
        throw std::runtime_error("Empty path in uninstall manifest.");

    const int count = MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS,
        relativeUtf8.data(), static_cast<int>(relativeUtf8.size()),
        nullptr, 0);

    if (count <= 0)
        throw std::runtime_error("Invalid UTF-8 path in uninstall manifest.");

    std::wstring text(static_cast<size_t>(count), L'\0');
    MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS,
        relativeUtf8.data(), static_cast<int>(relativeUtf8.size()),
        text.data(), count);

    std::replace(text.begin(), text.end(), L'\\', L'/');

    while (!text.empty() && text.front() == L'/')
        text.erase(text.begin());

    const fs::path raw(text);
    if (raw.is_absolute() || raw.has_root_name() ||
        raw.has_root_directory())
        throw std::runtime_error("Rooted path in uninstall manifest.");

    fs::path clean;
    for (const auto& part : raw) {
        const std::wstring component = part.wstring();

        if (component == L".." ||
            component.find(L':') != std::wstring::npos)
            throw std::runtime_error("Unsafe path in uninstall manifest.");

        if (!component.empty() && component != L".")
            clean /= part;
    }

    if (clean.empty())
        throw std::runtime_error("Empty path in uninstall manifest.");

    return root / clean;
}

int main() {
    try {
        const fs::path installDir = executable_path().parent_path();
        const fs::path manifest = installDir / L"uninstall-manifest.txt";
        const std::wstring appName = app_name();

        std::wcout << L"Uninstall McPhase from "
                   << installDir.wstring() << L"? (y/N) ";

        std::wstring answer;
        std::getline(std::wcin, answer);
        if (answer.empty() || (answer[0] != L'y' && answer[0] != L'Y')) {
            std::wcout << L"Uninstall cancelled.\n";
            return 0;
        }

        // Remove the shortcuts created by this installer.
        const fs::path programs =
            known_folder(FOLDERID_Programs) / appName;

        std::error_code ec;
        fs::remove(programs / (appName + L".lnk"), ec);
        fs::remove(programs / (appName + L" - Demo.lnk"), ec);
        fs::remove(programs / (appName + L" - Uninstall.lnk"), ec);
        fs::remove(programs, ec);

        const fs::path desktop = known_folder(FOLDERID_Desktop);
        fs::remove(desktop / (appName + L".lnk"), ec);
        fs::remove(desktop / (appName + L" - Demo.lnk"), ec);

        std::ifstream in(manifest, std::ios::binary);
        if (!in)
            throw std::runtime_error("Cannot open uninstall manifest.");

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty())
                continue;

            const fs::path target = safe_destination(installDir, line);
            fs::remove(target, ec);

            // Remove empty parent directories, stopping at the install folder.
            fs::path parent = target.parent_path();
            while (parent != installDir) {
                ec.clear();
                if (!fs::remove(parent, ec) || ec)
                    break;
                parent = parent.parent_path();
            }
        }

        fs::remove(manifest, ec);

        // Windows cannot delete the currently running executable.
        const fs::path self = executable_path();
        if (!MoveFileExW(self.c_str(), nullptr, MOVEFILE_DELAY_UNTIL_REBOOT))
            std::wcerr << L"Could not schedule uninstall.exe for removal.\n";

        std::wcout << L"Uninstall complete. Restart Windows to remove "
                      L"uninstall.exe.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Uninstall failed: " << e.what() << '\n';
        return 1;
    }
}
