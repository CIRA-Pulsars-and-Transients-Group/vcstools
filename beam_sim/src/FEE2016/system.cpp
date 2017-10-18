#include "system.h"

#include <boost/filesystem.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>

std::string System::FindPythonFilePath(const std::string& filename)
{
	if(boost::filesystem::exists(filename))
		return filename;
	std::cout << "Searching " << filename << "... " << std::flush;
	std::random_device rndDev;
  std::mt19937 gen(rndDev());
	std::stringstream filenameStr;
	filenameStr << "/tmp/mwa-ao-python-path-list" << gen() << ".tmp";
	boost::filesystem::path tempPath = filenameStr.str(); //boost::filesystem::unique_path();
	const std::string tempFilename = tempPath.string();  // optional
	std::string command =
		std::string("echo \"import sys\nfor a in sys.path:\n  print a\"|python>") +
		tempFilename;
	int status = system(command.c_str());
	if(status != 0)
	  throw std::runtime_error("system() returned non-zero error code: might be out of memory, or python might not be working properly");
	std::ifstream searchPathsFile(tempFilename.c_str());
	if(!searchPathsFile.good())
	  throw std::runtime_error(("Error in findPythonFilePath: system call did not create expected temporary file " + tempFilename).c_str());
	while(searchPathsFile.good())
	{
		std::string prefixPath;
		std::getline(searchPathsFile, prefixPath);
		boost::filesystem::path searchPath(prefixPath);
		searchPath /= filename;
		
		bool pathExists = false;
		try {
			pathExists = boost::filesystem::exists(searchPath);
		} catch(...) { }
		if(pathExists)
		{
			const std::string result = searchPath.string();
			std::cout << result << '\n';
			searchPathsFile.close();
			boost::filesystem::remove(tempPath);
			return result;
		}
	}
	searchPathsFile.close();
	boost::filesystem::remove(tempPath);
	throw std::runtime_error(std::string("Could not find Python file ") + filename);
}
