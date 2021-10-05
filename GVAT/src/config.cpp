#include "config.h"

using namespace std;

//constructor
Config::Config(string const& config_filename)
{
	ifstream infile(config_filename.c_str());	//open file
	if(!infile.is_open())	//open failed
		throw 99;
	
	string line;
	int line_num=1;
	while (getline(infile, line))	//lets read one line at a time then parse each line
	{
		istringstream ss(line);	//istringstream for parsing fields
		if(line[0] != '/')	//ignore lines that start with /
			ParseLine(ss);
		line_num++;
	}
	infile.close();
}

const string Config::GetFilename(string const& key)
{
	map<string,string>::iterator it;
	it = files_.find(key);
	if(it != files_.end())
		return it->second;
	else
		return "";
}

void Config::ParseLine(istringstream& ss)
//get a line
//parse the line and add the key and associated filename to our private map
{	
	string key, filename;
	ss >> key; //read key first
	ss >> filename; //then read path/filename
	files_.insert ( pair<string,string>(key, filename) );
	
}