#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <map>
#include <fstream>
#include <sstream>

using namespace std;
 
class Config{

public:
	Config(string const& config_filename);
	const string GetFilename(string const& key);
	
private:
	map <string,string> files_;
	void ParseLine(istringstream& ss);

};



#endif

