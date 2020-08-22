#pragma once
#include "lib.h"

class Client
{
public:
	vector<int> listLoc;
	int demand = 0;
	Client();
	~Client();

private:

};

Client::Client()
{
}

Client::~Client()
{
	listLoc.clear();
}