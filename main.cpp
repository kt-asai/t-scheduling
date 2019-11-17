#include <iostream>
#include <filesystem>

int main()
{
    std::cout << "This program is..." << std::endl;

    std::filesystem::path path = "benchmarks/tof_3.qc";

    std::cout << path.filename() << std::endl;

    return 0;
}
