#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct FileEntry {
    fs::path source;
    std::string destination; // Relative to {app}
};

static std::string trim(std::string s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// Split a directive at semicolons that are outside quoted strings.
static std::vector<std::string> split_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::string part;
    bool quoted = false;

    for (char c : line) {
        if (c == '"') quoted = !quoted;
        if (c == ';' && !quoted) {
            fields.push_back(trim(part));
            part.clear();
        } else {
            part += c;
        }
    }
    if (!trim(part).empty()) fields.push_back(trim(part));
    return fields;
}

static std::string unquote(std::string s) {
    s = trim(s);
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        return s.substr(1, s.size() - 2);
    return s;
}

static std::map<std::string, std::string> parse_fields(const std::string& line) {
    std::map<std::string, std::string> result;

    for (const auto& field : split_fields(line)) {
        const auto colon = field.find(':');
        if (colon == std::string::npos) continue;

        const auto key = lower(trim(field.substr(0, colon)));
        const auto value = unquote(field.substr(colon + 1));
        result[key] = value;
    }
    return result;
}

static bool has_flag(const std::string& flags, const std::string& wanted) {
    const std::string value = lower(flags);
    std::size_t pos = 0;

    while (pos < value.size()) {
        while (pos < value.size() &&
               (std::isspace(static_cast<unsigned char>(value[pos])) ||
                value[pos] == ',')) {
            ++pos;
        }

        const std::size_t start = pos;

        while (pos < value.size() &&
               !std::isspace(static_cast<unsigned char>(value[pos])) &&
               value[pos] != ',') {
            ++pos;
        }

        if (value.substr(start, pos - start) == wanted)
            return true;
    }

    return false;
}

static std::string normalize_separators(std::string s) {
    std::replace(s.begin(), s.end(), '\\', fs::path::preferred_separator);
    std::replace(s.begin(), s.end(), '/', fs::path::preferred_separator);
    return s;
}



static void add_file(std::vector<FileEntry>& files,
                     const fs::path& source,
                     const std::string& destination) {
    if (!fs::is_regular_file(source))
        throw std::runtime_error("Not a regular file: " + source.string());

    files.push_back({source, destination});
}

static bool wildcard_match(const std::string& pattern, const std::string& name) {
    // In Inno Setup-style file lists, *.* should include extensionless files.
    if (pattern == "*.*" || pattern == "*")
        return true;

    std::size_t p = 0, n = 0, star = std::string::npos, retry = 0;

    while (n < name.size()) {
        if (p < pattern.size() &&
            (pattern[p] == '?' ||
             std::tolower(static_cast<unsigned char>(pattern[p])) ==
                 std::tolower(static_cast<unsigned char>(name[n])))) {
            ++p;
            ++n;
        } else if (p < pattern.size() && pattern[p] == '*') {
            star = p++;
            retry = n;
        } else if (star != std::string::npos) {
            p = star + 1;
            n = ++retry;
        } else {
            return false;
        }
    }

    while (p < pattern.size() && pattern[p] == '*') ++p;
    return p == pattern.size();
}

static void collect_source(std::vector<FileEntry>& files,
                           const fs::path& project,
                           const std::string& source_text,
                           const std::string& dest_dir,
                           bool recurse) {
    const fs::path source_pattern(normalize_separators(source_text));
    const fs::path relative_parent = source_pattern.parent_path();
    const std::string filename_pattern = source_pattern.filename().string();

    const fs::path search_dir =
        (project / relative_parent).lexically_normal();

    if (!fs::exists(search_dir) || !fs::is_directory(search_dir))
        throw std::runtime_error("Source directory not found: " +
                                 search_dir.string());

    auto add_matching = [&](const fs::directory_entry& item) {
        if (!item.is_regular_file()) return;

        const std::string filename = item.path().filename().string();
        if (!wildcard_match(filename_pattern, filename)) return;

        // Preserve nested paths relative to this Source folder:
        // demo/sub/file.dat becomes {app}\demo\sub\file.dat.
        const fs::path relative_file =
            item.path().lexically_relative(search_dir);

        if (relative_file.empty() || relative_file.is_absolute())
            throw std::runtime_error("Could not make source path relative: " +
                                     item.path().string());

        std::string destination = dest_dir;
        if (!destination.empty() && destination.back() != '\\' &&
            destination.back() != '/') {
            destination += '\\';
        }

        destination += relative_file.generic_string();
        std::replace(destination.begin(), destination.end(), '/', '\\');

        add_file(files, item.path(), destination);
    };

    if (recurse) {
        // This iterator includes top-level files as well as files in
        // subdirectories, so don't also run directory_iterator here.
        for (const auto& item : fs::recursive_directory_iterator(search_dir))
            add_matching(item);
    } else {
        for (const auto& item : fs::directory_iterator(search_dir))
            add_matching(item);
    }
}


static std::vector<FileEntry> read_iss(const fs::path& script) {
    std::ifstream in(script);
    if (!in) throw std::runtime_error("Cannot open script: " + script.string());

    const fs::path project = script.parent_path();
    std::vector<FileEntry> files;
    std::string line;
    bool in_files_section = false;
    std::size_t line_number = 0;

    while (std::getline(in, line)) {
        ++line_number;
        const auto text = trim(line);

        if (text.empty() || text[0] == ';') continue;

        if (text.front() == '[' && text.back() == ']') {
            in_files_section = lower(text) == "[files]";
            continue;
        }
        if (!in_files_section) continue;

        try {
            const auto fields = parse_fields(text);
            const auto source = fields.find("source");
            const auto dest = fields.find("destdir");

            if (source == fields.end() || dest == fields.end()) {
                throw std::runtime_error(
                    "Each active [Files] entry needs Source and DestDir.");
            }

            std::string dest_dir = dest->second;
            const std::string app_prefix = "{app}";
            if (dest_dir.rfind(app_prefix, 0) == 0)
                dest_dir.erase(0, app_prefix.size());

            const auto flags = fields.find("flags");
            const bool recurse = flags != fields.end() &&
                                 has_flag(flags->second, "recursesubdirs");

            collect_source(files, project, source->second, dest_dir, recurse);
        } catch (const std::exception& e) {
            throw std::runtime_error("Line " + std::to_string(line_number) +
                                     ": " + e.what());
        }
    }

    if (files.empty())
        throw std::runtime_error("No files found in active [Files] entries.");

    return files;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: iss_files <script.iss> [output-manifest.tsv]\n";
        return 2;
    }

    
    try {
        const fs::path script = fs::absolute(argv[1]);
        const auto files = read_iss(script);

        const fs::path manifest =
        argc == 3 ? std::filesystem::path(argv[2])
                  : std::filesystem::path("files-manifest.tsv");

        std::ofstream out(manifest);
        if (!out) throw std::runtime_error("Cannot write: " + manifest.string());

        for (const auto& file : files) {
            out << file.source.generic_string() << '\t'
                << file.destination << '\n';
        }

        std::cout << "Collected " << files.size() << " files.\n"
                  << "Manifest: " << manifest << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
