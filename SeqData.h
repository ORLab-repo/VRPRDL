#pragma once

class SeqData
{
public:
	int firstnode;
	int lastnode;
	int cost;
	int load;
	int E;//earliest time completion
	int L;//latest starting time
	int T;//sum of travel and service time
	SeqData();
	~SeqData();

private:

};

SeqData::SeqData()
{
}

SeqData::~SeqData()
{
}