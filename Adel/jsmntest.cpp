// hg clone http://bitbucket.org/zserge/jsmn jsmn
// https://zserge.com/jsmn.html

#include <iostream>
#include "jsmn.h"
#include <cstring>

using namespace std;

jsmn_parser parser;
jsmntok_t* tokens;
char* js;

int main() {
    jsmn_init(&parser);
    js = new char[100];
    strncpy(js, "{\"name\":\"Ade\",\"age\": 20}", 30);
    int r = jsmn_parse(&parser, js, strlen(js), NULL, 0);
    if (r == JSMN_ERROR_INVAL) cout << "\033[1;31mbad token, JSON string is corrupted\033[0m" << endl;
    if (r == JSMN_ERROR_NOMEM) cout << "\033[1;31mnot enough tokens, JSON string is too large\033[0m" << endl;
    if (r == JSMN_ERROR_PART) cout << "\033[1;31mJSON string is too short, expecting more JSON data\033[0m" << endl;
    // tokens = new jsmntok_t[r];
    // jsmn_parse(&parser, js, strlen(js), tokens, r);
    return 0;
}