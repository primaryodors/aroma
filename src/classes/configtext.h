#ifndef AROMA_CONFIGTEXT_H
#define AROMA_CONFIGTEXT_H

#include <cstring>
#include <stdexcept>
#include <string>
#include <array>
#include <cstdlib>
#include <cmath>
#include <cctype>

inline std::array<float,3> parse_pdb_coordinates(const char* line)
{
    if (!line || std::strlen(line) < 54)
        throw std::invalid_argument("PDB coordinate record is too short.");
    std::array<float,3> result;
    for (std::size_t axis=0; axis<3; ++axis)
    {
        char field[9] = {};
        std::memcpy(field,line+30+8*axis,8);
        char* end = nullptr;
        result[axis] = std::strtof(field,&end);
        if (end == field || !std::isfinite(result[axis]))
            throw std::invalid_argument("Invalid PDB coordinate field.");
        while (*end && std::isspace(static_cast<unsigned char>(*end))) ++end;
        if (*end) throw std::invalid_argument("Invalid PDB coordinate field.");
    }
    return result;
}

template <std::size_t N>
void copy_config_text(char (&destination)[N], const char* source)
{
    if (!source || std::strlen(source) >= N)
        throw std::invalid_argument("Configuration text is missing or too long.");
    std::memcpy(destination, source, std::strlen(source)+1);
}

inline const char* config_basename(const char* filename)
{
    if (!filename) throw std::invalid_argument("Missing filename.");
    const char* slash = std::strrchr(filename, '/');
    const char* backslash = std::strrchr(filename, '\\');
    if (backslash && (!slash || backslash > slash)) slash = backslash;
    return slash ? slash+1 : filename;
}

inline std::string config_stem(const char* filename)
{
    std::string result=config_basename(filename);
    const auto dot=result.find('.');
    if(dot!=std::string::npos)result.resize(dot);
    return result;
}

#endif
