#pragma once
#include<vector>

using namespace std;

class Location
{
public:
	double x;
	double y;
	int stTime;//stating time
	int enTime;//ending time
	int idxClient = -1;
	vector<int> moves;

	Location()
	{
	}

	Location(double _x, double _y) {
		x = _x;
		y = _y;
	}

	int calDis(Location& val) {
		return round(sqrt((x - val.x) * (x - val.x) + (y - val.y) * (y - val.y)));
	}

	~Location()
	{
		moves.clear();
	}

private:
};
