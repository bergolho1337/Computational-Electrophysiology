#include "../include/parser.h"

Parser::Parser (int argc, char *argv[])
{
    string str;
    in.open(argv[1]);
    while (getline(in,str))
    {
        cout << str << endl;
        size_t found = str.find("ALGEBRAIC[");
        if (found != string::npos)
            cout << str << endl;
    } 
    in.close();
}

void Parser::convert ()
{
    cout << "[+] Converting ..." << endl;
}