#pragma once
#include "lib.h"

class Client
{
public:
	vector<int> listLoc;
	int demand = -1;
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