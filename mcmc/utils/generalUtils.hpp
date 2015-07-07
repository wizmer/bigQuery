#ifndef _GENERALUTILS__HPP_
#define _GENERALUTILS__HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <map>
#include <vector>
#include <algorithm>

#include <stdio.h>

#include "dirent.h"

namespace generalUtils{
    std::map< std::string, float > getEnergyBeamTest();
    std::map< std::string, std::string > getEnergyBeamTestString();
    int getColor(int n, bool fill = 0);

    std::string getPath(std::string fullFileName);
    std::string getFileName(std::string fullFileName);
    std::string getExtension( std::string fileName );

    std::vector <std::string > getFilesInDir(std::string dirName, int maxNumberOfFile = 0);
    std::vector <std::string > getFilesInDirWithPattern(std::string dirName, std::string pattern);

    std::string dateTime();
    bool folderExists( std::string folderName );
    bool fileExists( std::string fileName );

    int makeFolder( std::string folderName );

    // return STDOUT from the shell command: cmd (but not STDERR)
    std::string exec(std::string cmd);

    // split a string containing line returns into several strings
    std::vector<std::string> splitIntoLines(const std::string &string);
    std::vector< std::string > split( std::string str, std::string delimiters );
    // Same as Split but with reverse ordering
    std::vector< std::string > splitFromEnd( std::string str, std::string delimiters );
    std::string replacePattern(std::string str, const std::string& oldPattern, const std::string& newPattern);

    template<class T> T stringTo(std::string str){
        std::stringstream ss;
        T res;
        ss << str;
        ss >> res;
        return res;
    }

    template<class T> std::vector<T> stringTo(std::vector<std::string> stringVector){
        int size = stringVector.size();
        std::vector<T> output;
        for(int i = 0;i<size;i++){
            std::stringstream ss;
            T res;
            ss << stringVector[i];
            ss >> res;
            output.push_back(res);
	}
        return output;
    }

    std::string toString(float a);
}

#endif

