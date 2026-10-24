#ifndef AROMA_CONFIGTEXT_H
#define AROMA_CONFIGTEXT_H

#include <cstring>
#include <stdexcept>
#include <string>

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

#endif
