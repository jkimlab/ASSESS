#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>
#include <list>
#include <iomanip>
#include <cmath>

using namespace std;

enum {SUCCESS, FAIL, CYCLE}; 
double gMinAdjScore = 0;

class Edge {
public:
	int bid1, bid2;
	int dir1, dir2;
	double score1, score2, weight;

	Edge(int _bid1, int _dir1, int _bid2, int _dir2)
			:bid1(_bid1), dir1(_dir1), bid2(_bid2), dir2(_dir2), 
			score1(0.0), score2(0.0), weight(0.0) {}
	
	bool operator<(const Edge& p) const {
		if (bid1 == p.bid1 && dir1 == p.dir1 && bid2 == p.bid2) return dir2 < p.dir2;
		if (bid1 == p.bid1 && dir1 == p.dir1) return bid2 < p.bid2;
		if (bid1 == p.bid1) return dir1 < p.dir1;
		return bid1 < p.bid1;
	}

	void reverse() {
		int tmp = bid1;
		bid1 = bid2;
		bid2 = tmp;

		tmp = dir1;
		dir1 = dir2;
		dir2 = tmp;
		if (dir1 == 1) dir1 = -1;
		else dir1 = 1;
		if (dir2 == 1) dir2 = -1;
		else dir2 = 1;
	}

	string toString() {
		stringstream ss;
		ss << bid1 * dir1 << " " << bid2 * dir2;
		return ss.str();
	}
};

void error (string msg, string file="") 
{
	cerr << msg << file << endl;
	exit(1);
}

bool cmp(const pair<Edge,double> &p1, const pair<Edge,double> &p2) 
{
    if (p1.second > p2.second) {
        return true;
    } else if (p1.second == p2.second) {
        int dist1 = abs(p1.first.bid1 - p1.first.bid2); 
        int dist2 = abs(p2.first.bid1 - p2.first.bid2);
        if (dist1 < dist2) return true;
        else return false;
    } else {
        return false;
    }
    //return p1.second > p2.second;
}

int insertEdge(list<Edge>& le, Edge& e, map<int,int>& mapUsed)
{
	Edge& fe = le.front();
	Edge& be = le.back();
	
	if (fe.bid1 == e.bid1 && fe.dir1 != e.dir1) {
		// check for a cycle
		if (be.bid2 == e.bid2 && be.dir2 != e.dir2) {
			return CYCLE; 
		}
		
		// e precedes fe with an opposite direction
		e.reverse();
		le.push_front(e);
		mapUsed[fe.bid1] = 1;
		mapUsed[-fe.bid1] = 1;
		if (e.dir1 == 1) mapUsed[-e.bid1] = 1;
		else mapUsed[e.bid1] = 1;
		return SUCCESS;	
	}
	if (fe.bid1 == e.bid2 && fe.dir1 == e.dir2) {
		
		// check for a cycle
		if (be.bid2 == e.bid1 && be.dir2 == e.dir1) return CYCLE; 

		// e precedes fe with the same direction
		le.push_front(e);
		mapUsed[fe.bid1] = 1;
		mapUsed[-fe.bid1] = 1;
		if (e.dir1 == 1) mapUsed[-e.bid1] = 1;
		else mapUsed[e.bid1] = 1;
		return SUCCESS;
	}

	if (be.bid2 == e.bid1 && be.dir2 == e.dir1) {
		// check for a cycle
		if (fe.bid1 == e.bid2 && fe.dir1 == e.dir2) return CYCLE; 
		
		// be precedes e with the same direction
		le.push_back(e);
		mapUsed[be.bid2] = 1;
		mapUsed[-be.bid2] = 1;
		if (e.dir2 == 1) mapUsed[e.bid2] = 1;
		else mapUsed[-e.bid2] = 1;
		return SUCCESS;
	}
	if (be.bid2 == e.bid2 && be.dir2 != e.dir2) {
		// check for a cycle
		if (fe.bid1 == e.bid1 && fe.dir1 != e.dir1) return CYCLE; 
		
		// be precedes e with an opposite direction
		e.reverse();
		le.push_back(e);
		mapUsed[be.bid2] = 1;
		mapUsed[-be.bid2] = 1;
		if (e.dir2 == 1) mapUsed[e.bid2] = 1;
		else mapUsed[-e.bid2] = 1;
		return SUCCESS;	
	}
	return FAIL;
}

void printLists(map<int, list<Edge> >& mapClasses)
{
	map<int, list<Edge> >::iterator citer;

	int clsnum = 1;	
	for(citer = mapClasses.begin(); citer != mapClasses.end(); citer++) {
		list<Edge>& le = citer->second;
		cout << "#" << clsnum << endl;
		clsnum++;
	
		list<Edge>::iterator liter;
		for (liter = le.begin(); liter != le.end(); liter++) {
			Edge& e = *liter;
            
            if (next(liter) == le.end()) {
			    cout << e.bid1 * e.dir1 << " " << e.bid2 * e.dir2 << " $"; 
            } else {
			    cout << e.bid1 * e.dir1 << " "; 
            }
		}	
		cout << endl;
	} // end of for
}

void mergeLists(int clscnt, int clsid, list<Edge>& le, map<int, list<Edge> > &mapClasses) 
{
	map<int, list<Edge> >::iterator citer;
	list<Edge>::iterator liter;
	list<Edge>::reverse_iterator rliter;
	list<Edge>& le1 = le;
		
	for (int j = 1; j <= clscnt; j++) {
		if (j == clsid) continue;

		citer = mapClasses.find(j);
		if (citer == mapClasses.end()) continue;
		list<Edge>& le2 = citer->second;
			
		Edge& e1front = le1.front();
		Edge& e1back = le1.back();
	
		Edge& e2front = le2.front();
		Edge& e2back = le2.back();

		if(e1front.bid1 == e2front.bid1 && e1front.dir1 != e2front.dir1) {
			// check cycle
			if (e1back.bid2 != e2back.bid2) {	
				for (liter = le2.begin(); liter != le2.end(); liter++) {
					Edge e2 = *liter;
					e2.reverse();
					le1.push_front(e2);
				} // end of for
				mapClasses.erase(j);
			}
		} else if(e1front.bid1 == e2back.bid2 && e1front.dir1 == e2back.dir2) {
			// check cycle
			if (e1back.bid2 != e2front.bid1) {
				for (rliter = le2.rbegin(); rliter != le2.rend(); rliter++) {
					Edge e2 = *rliter;
					le1.push_front(e2);
				} // end of for
				mapClasses.erase(j);	
			}	
		} else if(e1back.bid2 == e2front.bid1 && e1back.dir2 == e2front.dir1) {
			// check cycle
			if (e1front.bid1 != e2back.bid2) {
				for (liter = le2.begin(); liter != le2.end(); liter++) {
					Edge e2 = *liter;
					le1.push_back(e2);
				} // end of for
				mapClasses.erase(j);	
			}	
		} else if(e1back.bid2 == e2back.bid2 && e1back.dir2 != e2back.dir2) {
			// check cycle
			if (e1front.bid1 != e2front.bid1) {
				for (rliter = le2.rbegin(); rliter != le2.rend(); rliter++) {
					Edge e2 = *rliter;
					e2.reverse();
					le1.push_back(e2);
				} // end of for
				mapClasses.erase(j);	
			}
		}	
	} // end of for j
}

int main(int argc, char* argv[]) 
{
	string fjoinbreak = string(argv[1]);
	gMinAdjScore = strtod(argv[2], NULL);
	char* fscore = argv[3];

    cerr << "Minimum adjacency score = " << gMinAdjScore << endl;	
	cerr << "Forced join/break file = " << fjoinbreak << endl;

	// read forced join/break data
	map<Edge, string> mapFJoinBreak;
	ifstream infile;
	if (fjoinbreak != "none") {
		infile.open(fjoinbreak.c_str());
		if (!infile) error ("\n[ERROR] Unable to open file: ", fjoinbreak);
		int bid1, bid2;
		string label;
		while (infile >> bid1 >> bid2 >> label) {
			int intdir1 = (bid1 > 0) ? 1 : -1;
			int intdir2 = (bid2 > 0) ? 1 : -1;
			Edge bpf(abs(bid1), intdir1, abs(bid2), intdir2);
			Edge bpr(abs(bid2), -intdir2, abs(bid1), -intdir1);
			mapFJoinBreak[bpf] = label;
			mapFJoinBreak[bpr] = label;
		} // end of while
		infile.close();	
	}

	// read adjacency scores
    map<Edge, string>::iterator esiter;
	map<Edge, double> mapAdjScores;
	infile.open(fscore);
	if (!infile) error ("\n[ERROR] Unable to open file: ", fscore); 
	int bid1, bid2;
	double adjscore;
	while (infile >> bid1 >> bid2 >> adjscore) {
        int intdir1 = (bid1 > 0) ? 1 : -1;
        int intdir2 = (bid2 > 0) ? 1 : -1;
        Edge bpf(abs(bid1), intdir1, abs(bid2), intdir2);
        Edge bpr(abs(bid2), -intdir2, abs(bid1), -intdir1);

        esiter = mapFJoinBreak.find(Edge(abs(bid1), intdir1, abs(bid2), intdir2));
        if (esiter != mapFJoinBreak.end()) {
            string label = esiter->second;
            if (label == "JOIN") {
                adjscore = 2;
            } else if (label == "BREAK") {
                adjscore = 0;
            } else {
                error("\n[ERROR] Invalid label in the forced join/break file: ", label);
            }

            mapFJoinBreak.erase(esiter);
        }
        
        if (adjscore == 0) continue;
		
        mapAdjScores[bpf] = adjscore;
		mapAdjScores[bpr] = adjscore;
	}
	infile.close();

    for (esiter = mapFJoinBreak.begin(); esiter != mapFJoinBreak.end(); esiter++) {
        Edge bpf = esiter->first;
        string label = esiter->second;
        double adjscore;  

        if (label == "JOIN") {
            adjscore = 2;
        } else if (label == "BREAK") {
            adjscore = 0;
        } else {
            error("\n[ERROR] Invalid label in the forced join/break file: ", label);
        }
        
        Edge bpr(bpf.bid2, -1 * bpf.dir2, bpf.bid1, -1 * bpf.dir1);
        
        mapAdjScores[bpf] = adjscore;
        mapAdjScores[bpr] = adjscore;
    }

	// greedy search based on edge weights
	vector<pair<Edge,double> > vecEdges(mapAdjScores.begin(), mapAdjScores.end()); 
	sort (vecEdges.begin(), vecEdges.end(), cmp);

	map<int, int> mapUsed;	// pos:head is connected, neg:tail is connected
	map<int, int>::iterator uiter, uiterex;
	map<int, list<Edge> > mapClasses;
	map<int, list<Edge> >::iterator citer;
	map<Edge, int>::iterator niter;
	map<Edge, double>::iterator biter;
	int clscnt = 0;
	for (int i = 0; i < vecEdges.size(); i++) {
		pair<Edge,double> p = vecEdges.at(i);
		Edge& e = p.first;
		e.weight = p.second;

		if (e.weight < gMinAdjScore) continue;

		if (e.dir1 == 1) uiter = mapUsed.find(-e.bid1); 
		else uiter = mapUsed.find(e.bid1);  
		if (uiter != mapUsed.end()) {
			continue;				
		}
	
		if (e.dir2 == 1) uiter = mapUsed.find(e.bid2);
		else uiter = mapUsed.find(-e.bid2); 
		if (uiter != mapUsed.end()) {
			continue;		
		}		
		
        bool found = false;
		for (citer = mapClasses.begin(); citer != mapClasses.end(); citer++) {
			list<Edge>& le = citer->second;
			int res = insertEdge(le, e, mapUsed);
			if (res == SUCCESS || res == CYCLE) { 
				if (res == SUCCESS) mergeLists(clscnt, citer->first, le, mapClasses);
				found = true;
				break;
			} // end of if	
		} // end of for 
	
		if (found == false) {
			list<Edge> le;
			le.push_back(e);
			mapClasses[++clscnt] = le;
		
			if (e.dir1 == 1) mapUsed[-e.bid1] = 1; 
			else mapUsed[e.bid1] = 1; 
			if (e.dir2 == 1) mapUsed[e.bid2] = 1; 
			else mapUsed[-e.bid2] = 1; 
		}
	}

	printLists(mapClasses);

	return 0;
}
