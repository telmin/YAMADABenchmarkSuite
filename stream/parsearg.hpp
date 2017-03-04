#ifndef PARSE_ARG_HPP
#define PARSE_ARG_HPP

#include <iostream>
#include <string>
#include <cstring>

class ParseArg
{
public:
    ParseArg(int argc, char** argv)
	: m_array_size(100000), m_times(10)
    {
	parse(argc, argv);
    }

    size_t getArraySize() const {
	return m_array_size;
    }

    size_t getTimes() const {
	return m_times;
    }

private:
    bool convert(const std::string& str, size_t& ret) {
	char* ll;
	unsigned int r = strtoul(str.c_str(), &ll, 10);
	ret = r;
	return !strlen(ll);
    }

    bool convertSuffix(const std::string& str, size_t& ret) {
	char* ll;
	unsigned int r = strtoul(str.c_str(), &ll, 10);

	unsigned int coeff = 1;
	if(strlen(ll) != 0) {
	    switch(ll[0]) {
	    case 'G':
	    case 'g':
		coeff = 1024 * 1024 * 1024;
		break;
	    case 'M':
	    case 'm':
		coeff = 1024 * 1024;
		break;
	    case 'K':
	    case 'k':
		coeff = 1024;
		break;
	    default:
		break;
	    }
	}

	ret = static_cast<size_t>(r) * coeff;

	return true;
    }


    void parse(int argc, char** argv) {
	for(int i = 1; i < argc; ++i) {
	    if(!std::string("--numtimes").compare(argv[i]) ||
	       !std::string("-n").compare(argv[i])) {
		if(++i >= argc || !convert(argv[i], m_times)) {
		    std::cerr << "Invalid number of times." << std::endl;
		    exit(EXIT_FAILURE);
		}
		if(m_times < 2) {
		    std::cout << "!!warning!! ntimes must be greater than 1." << std::endl;
		    m_times = 10;
		}
	    } else if(!std::string("--arraysize").compare(argv[i]) ||
		      !std::string("-s").compare(argv[i])) {
		if(++i >= argc || !convertSuffix(argv[i], m_array_size)) {
		    std::cerr << "Invalid array size." << std::endl;
		    exit(EXIT_FAILURE);
		}
	    } else if(!std::string("--byte").compare(argv[i]) ||
		      !std::string("-b").compare(argv[i])) {
		size_t bytesize;
		if(++i >= argc || !convertSuffix(argv[i], bytesize)) {
		    std::cerr << "Invalid byte size." << std::endl;
		    exit(EXIT_FAILURE);
		}

		m_array_size = bytesize / sizeof(double);
	    }
	}
    }

    size_t m_array_size;
    size_t m_times;
};


#endif
