#pragma once

class SeqData
{
public:
	int firstnode;//index location
	int lastnode;//index location
	int cost;
	int load;
	int E;//earliest time completion
	int L;//latest starting time
	int T;//sum of travel and service time
	bool F;//check feasibility;
	SeqData() {};
	~SeqData() {};

private:

};