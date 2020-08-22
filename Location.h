#pragma once
#include<vector>

using namespace std;

class Location
{
public:
	int x;
	int y;
	int stTime;//stating time
	int enTime;//ending time
	int idxClient = -1;
	vector<int> moves;
	Location();
	Location(int _x, int _y);
	~Location();

private:

};

Location::Location()
{
}

Location::Location(int _x, int _y) {
	x = _x;
	y = _y;
}

Location::~Location()
{
	moves.clear();
}
