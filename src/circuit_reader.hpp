#ifndef T_SCHEDULING_CIRCUIT_READER_HPP
#define T_SCHEDULING_CIRCUIT_READER_HPP

#include <vector>
#include <sstream>

#include "circuit.hpp"


class CircuitReader
{
private:
    std::string path_;

    /**
     * extract extension from file path
     * @param path file path
     * @return extension
     */
    inline std::string getExtension_(const std::string& path)
    {
        std::string ext;
        size_t pos1 = path.rfind('.');
        if (pos1 != std::string::npos)
        {
            ext = path.substr(pos1+1, path.size()-pos1);
            std::string::iterator itr = ext.begin();
            while (itr != ext.end())
            {
                *itr = tolower(*itr);
                itr++;
            }
            itr = ext.end()-1;
            while (itr != ext.begin())
            {
                if (*itr == 0 || *itr == 32)
                {
                    ext.erase(itr--);
                }
                else
                {
                    itr--;
                }
            }
        }
        return ext;
    }

    /**
     * split string by delimiter
     * @param input a input string
     * @param delimiter a delimiter as like space and ',' etc...
     * @return some splitted strings
     */
    inline std::vector<std::string> splitString_(const std::string& input,
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
    void readQC_(Circuit& qc);

public:
    /**
     * constructor
     */
    CircuitReader() = default;

    /**
     * constructor
     * @param path file path
     */
    CircuitReader(const std::string path) : path_(path) { }

    Circuit read();
};

#endif //T_SCHEDULING_CIRCUIT_READER_HPP
