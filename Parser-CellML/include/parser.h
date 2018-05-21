#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class Element
{
private:
    string name;
    int id;
public:
    Element ();
};

class Parser
{
private:
    ifstream in;
    Element *elem;
public:
    Parser (int argc, char *argv[]);
    void convert ();
};