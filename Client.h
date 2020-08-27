#pragma once
#include<vector>

using namespace std;

class Client
{
public:
	vector<int> listLoc;
	int demand = -1;

	Client()
	{
	}

	~Client()
	{
		listLoc.clear();
	}
};
