#ifndef T_SCHEDULING_CIRCUIT_READER_HPP
#define T_SCHEDULING_CIRCUIT_READER_HPP

#include <string>
#include <filesystem>

class CircuitReader
{
private:
    std::string file_name_;



public:
    CircuitReader() = default;

    CircuitReader(const std::string file_name) : file_name_(file_name) { }

    Circuit read();

};

#endif //T_SCHEDULING_CIRCUIT_READER_HPP
