#ifndef T_SCHEDULING_CIRCUIT_READER_HPP
#define T_SCHEDULING_CIRCUIT_READER_HPP

#include <vector>

#include "circuit.hpp"


class CircuitReader
{
private:
    std::string path_;
    std::vector<int> qubit_names_;

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
            {    // パスの最後に\0やスペースがあったときの対策
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

    void readQC_(Circuit& qc);

public:
    CircuitReader() = default;

    CircuitReader(std::string path) : path_(path) { }

    Circuit read();
};

#endif //T_SCHEDULING_CIRCUIT_READER_HPP
