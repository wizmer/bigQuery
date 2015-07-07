#include "generalUtils.hpp"

namespace generalUtils{
    std::string redFontBegin = "\033[1;31m";
    std::string redFontEnd = "\033[0m";


    std::map< std::string, float > getEnergyBeamTest(){
        std::map< std::string, float > energyBeamTest;
        energyBeamTest["pos10"]  = 10;
        energyBeamTest["pos20"]  = 20;
        energyBeamTest["pos80"]  = 80;
        energyBeamTest["pos100"] = 100;
        energyBeamTest["pos120"] = 119.9;    
        energyBeamTest["pos180"] = 179.5;
        energyBeamTest["pos300"] = 300;

        energyBeamTest["el100"] = 99.9;
        energyBeamTest["el120"] = 119.9;
        energyBeamTest["el180"] = 178.6;
        energyBeamTest["el300"] = 290.8;

        return energyBeamTest;
    }

    std::map< std::string, std::string > getEnergyBeamTestString(){
        std::map< std::string, std::string > energyBeamTestString;
        std::map< std::string, float > energyBeamTest = getEnergyBeamTest();

        for( std::map< std::string, float>::iterator it = energyBeamTest.begin(); it != energyBeamTest.end(); ++it){
            std::stringstream ss;
            std::string energyString;
            ss << it -> second;
            ss >> energyString;
            energyBeamTestString[ it -> first ] = energyString;
        }

        return energyBeamTestString;
    }
  
    int getColor(int n, bool fill ){
        int col[8]={4, 2, 8, 1, 7, 46, 6, 49};
        int color=n<8?col[n]:5+n;
        if (n==0 && fill) color=38;
        return color;
    }

    std::string getPath(std::string fullFileName){
        //Get old file, old tree and set top branch address
        std::string path;
        std::string fileName;

        size_t pos = fullFileName.find_last_of('/');
        if(pos != std::string::npos){
            path = fullFileName.substr(0,pos+1);
        }
        return path;
    }  

    std::string getFileName(std::string fullFileName){
        std::string path;
        std::string fileName;

        size_t pos = fullFileName.find_last_of('/');
        if(pos == std::string::npos){
            fileName = fullFileName;
        }else{
            fileName = fullFileName.substr(pos+1, fullFileName.length());
        }
        return fileName;
    }

    std::string getExtension( std::string fileName ){
        std::string extension;
    
        size_t pos = fileName.find_last_of('.');
        if(pos != std::string::npos){
            extension = fileName.substr(pos+1, fileName.length());
        }
        return extension;
    }

    std::vector <std::string > getFilesInDir(std::string dirName, int maxNumberOfFile ){
        std::vector < std::string > files;
        DIR *dir;
        struct dirent *ent;
        int i = 0;
        if ((dir = opendir (dirName.c_str())) != NULL) {
            /* print all the files and directories within directory */
            while ((ent = readdir (dir)) != NULL) {
                if(maxNumberOfFile > 0 && i >= maxNumberOfFile) break;
                std::string file(ent->d_name);
                if(file == "." || file == "..") continue;
                files.push_back( dirName+'/'+file );
                i++;
            }
            closedir (dir);
        }
        return files;
    }

    std::vector <std::string > getFilesInDirWithPattern(std::string dirName, std::string pattern){
        std::vector < std::string > files = getFilesInDir( dirName );
        std::vector < std::string > filesWithPattern;
        for(int i = 0; i < files.size(); ++i){
            size_t found = files[i].find( pattern );
            if( found != std::string::npos ) filesWithPattern.push_back( files[i] );
        }
        return filesWithPattern;
    }

    std::string dateTime(){
        time_t now = time(0);
        tm *ltm = localtime(&now);

        std::stringstream ss;
        ss << 1900 + ltm->tm_year;
        if(1 + ltm->tm_mon < 10) ss << 0;
        ss << 1 + ltm->tm_mon;
        if( ltm->tm_mday < 10) ss << 0;
        ss << ltm->tm_mday << '_';
        if( ltm->tm_hour < 10) ss << 0;
        ss << ltm->tm_hour;
        if( ltm->tm_min < 10) ss << 0;
        ss <<  ltm->tm_min;
        if( ltm->tm_sec < 10) ss << 0;
        ss << ltm->tm_sec;
        return ss.str();
    }

    bool folderExists( std::string folderName ){
        struct stat st;
        stat(folderName.c_str() , &st);
        return S_ISDIR(st.st_mode);
    }

    bool fileExists( std::string fileName ){
        return( std::ifstream(fileName.c_str()) != NULL );
    }

    int makeFolder( std::string folderName ){
        return mkdir( folderName.c_str() , S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    // return STDOUT from the shell command: cmd (but not STDERR)
    std::string exec(std::string cmd) {
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) return "ERROR";
        char buffer[128];
        std::string result = "";
        while(!feof(pipe)) {
            if(fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
        pclose(pipe);
        return result;
    }

    std::vector<std::string> splitIntoLines(const std::string &string){
        std::stringstream stream(string);
        std::vector<std::string> res;
        while (1){
            std::string line;
            std::getline(stream,line);
            if (!stream.good())
                break;
            res.push_back(line);
        }
        return res;
    }

    std::vector< std::string > split( std::string str, std::string delimiters ){
        std::vector < std::string > res;
        if( str.length() == 0 ) return res;

        std::size_t pos;

        do{
            pos = str.find_first_of( delimiters );
            //    cout << "avant : " << str << endl; 
            std::string a = str.substr(0,pos);
            if( a.length() > 0) res.push_back(a);
            str = str.substr(pos+1);
            //    cout << "apres : " << str << endl;
        } while( pos != std::string::npos );

        return res;
    }

    // Split in reverse order
    std::vector< std::string > splitFromEnd( std::string str, std::string delimiters ){
        std::vector< std::string > splittedString = split( str, delimiters );
        std::reverse( splittedString.begin(), splittedString.end() );
        return splittedString;
    }

    std::string replacePattern(std::string str, const std::string& oldPattern, const std::string& newPattern)
    {
        std::string result = str;
        size_t pos = 0;
        while((pos = result.find(oldPattern, pos)) != std::string::npos)
            {
                result.replace(pos, oldPattern.length(), newPattern);
                pos += newPattern.length();
            }

        return result;
    }

    std::string toString(float a){
        std::stringstream ss;
        std::string res;
        ss << a;
        ss >> res;
        return res;
    }
}
