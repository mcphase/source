#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "miniz.h"

namespace fs = std::filesystem;

static void write_u32(std::ofstream& out, std::uint32_t value) {
    for (int i = 0; i < 4; ++i)
        out.put(static_cast<char>((value >> (i * 8)) & 0xff));
    if (!out) throw std::runtime_error("Failed writing archive.");
}

static void write_u64(std::ofstream& out, std::uint64_t value) {
    for (int i = 0; i < 8; ++i)
        out.put(static_cast<char>((value >> (i * 8)) & 0xff));
    if (!out) throw std::runtime_error("Failed writing archive.");
}

static std::vector<unsigned char> read_file(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read file: " + path.string());

    in.seekg(0, std::ios::end);
    const auto end = in.tellg();
    if (end < 0) throw std::runtime_error("Cannot determine file size: " + path.string());

    std::vector<unsigned char> data(static_cast<std::size_t>(end));
    in.seekg(0, std::ios::beg);

    if (!data.empty() &&
        !in.read(reinterpret_cast<char*>(data.data()),
                 static_cast<std::streamsize>(data.size()))) {
        throw std::runtime_error("Failed reading file: " + path.string());
    }

    return data;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: package <files-manifest.tsv> <payload.mcp>\n";
        return 2;
    }

    try {
        std::ifstream manifest(argv[1]);
        if (!manifest) throw std::runtime_error("Cannot open manifest.");

        struct Entry {
            fs::path source;
            std::string destination;
        };

        std::vector<Entry> entries;
        std::string line;

        while (std::getline(manifest, line)) {
            if (line.empty()) continue;

            const auto tab = line.find('\t');
            if (tab == std::string::npos)
                throw std::runtime_error("Invalid manifest line; expected a tab.");

            entries.push_back({
                fs::u8path(line.substr(0, tab)),
                line.substr(tab + 1)
            });
        }

        if (entries.empty()) throw std::runtime_error("Manifest contains no files.");
        if (entries.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("Too many files for this archive format.");

        std::ofstream out(argv[2], std::ios::binary);
        if (!out) throw std::runtime_error("Cannot create output archive.");

        out.write("MCP1", 4);
        write_u32(out, static_cast<std::uint32_t>(entries.size()));

        for (const auto& entry : entries) {
            const auto path_bytes = entry.destination;
            if (path_bytes.size() > std::numeric_limits<std::uint32_t>::max())
                throw std::runtime_error("Destination path is too long.");

            const auto source = read_file(entry.source);
            if (source.size() > std::numeric_limits<mz_ulong>::max())
                throw std::runtime_error("File is too large for this miniz build.");

            const mz_ulong source_size = static_cast<mz_ulong>(source.size());
            const mz_ulong capacity = mz_compressBound(source_size);
            std::vector<unsigned char> compressed(static_cast<std::size_t>(capacity));
            mz_ulong compressed_size = capacity;

            const int result = mz_compress2(
    compressed.data(),
    &compressed_size,
    source.empty() ? nullptr : source.data(),
    source_size,
    MZ_BEST_COMPRESSION
);


            if (result != MZ_OK)
                throw std::runtime_error("Compression failed for: " + entry.source.string());

            write_u32(out, static_cast<std::uint32_t>(path_bytes.size()));
            write_u64(out, static_cast<std::uint64_t>(source.size()));
            write_u64(out, static_cast<std::uint64_t>(compressed_size));

            out.write(path_bytes.data(),
                      static_cast<std::streamsize>(path_bytes.size()));
            out.write(reinterpret_cast<const char*>(compressed.data()),
                      static_cast<std::streamsize>(compressed_size));

            if (!out) throw std::runtime_error("Failed writing archive data.");

            std::cout << entry.source << " -> " << entry.destination
                      << " (" << source.size() << " -> "
                      << compressed_size << " bytes)\n";
        }

        std::cout << "Wrote " << entries.size() << " files to " << argv[2] << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}

