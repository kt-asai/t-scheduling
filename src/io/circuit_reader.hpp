#ifndef T_SCHEDULING_CIRCUIT_READER_HPP
#define T_SCHEDULING_CIRCUIT_READER_HPP

#include <vector>
#include <sstream>

#include "../circuit/circuit.hpp"

namespace tskd {

class CircuitReader
{
private:
    std::string path_;

    /**
     * extract extension from file path
     * @param path file path
     * @return extension
     */
    static inline std::string extension(const std::string& path)
    {
        std::string extension;
        std::string::size_type pos = path.find_last_of('.');
        if (pos != std::string::npos)
        {
            extension = path.substr(pos + 1, path.size());
        }
        return extension;
    }

    /**
     * split string by delimiter
     * @param input a input string
     * @param delimiter a delimiter as like space and ',' etc...
     * @return some splitted strings
     */
    static inline std::vector<std::string> split(const std::string& input,
                                                 char delimiter)
    {
        std::istringstream stream(input);
        std::string field;
        std::vector<std::string> result;
        while (getline(stream, field, delimiter))
        {
            result.push_back(field);
        }
        return result;
    }

    /**
     * read qc format
     * @param qc empty Circuit class
     */
    void read_qc(Circuit& circuit);

public:
    /**
     * constructor
     */
    CircuitReader() = default;

    /**
     * constructor
     * @param path file path
     */
    CircuitReader(const std::string& path)
        : path_(path) { }

    Circuit read();
};

}

#endif //T_SCHEDULING_CIRCUIT_READER_HPP
