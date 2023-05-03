#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <regex>

using namespace std;

#define n_visit 1
#define n_unvisit 0

vector<string> Split_by_space (const string& str);
vector<string> Return_reverse (string crd_a, string crd_b);

class _Node {
	bool visit;
	bool used;
	map<string, float> link;
  public:
	_Node(string id, float score) { visit = false; used = false; link[id] = score; };
	void Update_used(){ used = true; }
	void Update_visit(int state) { if(state){ visit = true; } else { visit = false; }}
	bool Is_visit () { return visit; }
	bool Is_used () { return used; }
	void Add_link (string id, float score);
	void Delete_link (string id) { link.erase(id); }

	map<string, float> Get_adjNode(){ return link; }
	void Print_link ();
};

inline void _Node::Add_link (string id, float score) { 
	if(this->link.find(id) == link.end()){
		this->link[id] = score; 
	}else{
		this->link[id] += score; 
	}
}

inline void _Node::Print_link (){
	for(map<string, float>::iterator it = link.begin(); it != link.end(); ++it){
		cerr  << "-> " << it->first << "(" << it->second << ")";
	}
	cerr << endl;
}

struct _Link {
	string src, tar;
	float score;
};

bool operator < (const _Link &l1, const _Link &l2){
	return l1.score > l2.score;
}

bool operator == (const _Link &l1, const _Link &l2){
	return (l1.src == l2.src && l1.tar == l2.tar);
}

class _Network {
	map <string, _Node*> node_list;
	vector <string> next_visit;

	friend class _Traveler;
  public:
	int Save_info (string src, string tar, float score);
	int Build_network (vector <string> files);
	bool Search_node (string node) { return (node_list.find(node) == node_list.end())? false : true; }
	int Search_localgreedypath (const string start, const string goal, float score, map<vector<string>, float>& path);
	map<string, vector<string> > Search_globalgreedypath (const vector<string>& bs, const vector<string>& es);

	void Print_network();
};

inline void _Network::Print_network(){
	for(map <string, _Node*>::iterator it = node_list.begin(); it != node_list.end(); ++it){
		map <string, float> tars = it->second->Get_adjNode();
		for(map <string, float>::iterator mit = tars.begin(); mit != tars.end(); ++mit){
			cout << it->first << "\t" << mit->first << "\t" << mit->second << endl;	
		}
	}
};

int _Network::Save_info (string src, string tar, float score){
	if(node_list.find(src) == node_list.end()){
		node_list.insert(make_pair(src, new _Node(tar, score)));			
	}else{
		node_list[src]->Add_link(tar, score);
	}
}

int _Network::Build_network (vector <string> files){
	string ln;

	int file_cnt = files.size();

	for(int i=0; i<file_cnt; ++i){
		string ln, file = files[i];
		cerr << "Read " << file << " file" << endl;
		
		ifstream openFile(file.data());
		if(!openFile.is_open()){
			cerr << "[Error] Cannot open " << file << endl;
			return 0;
		}

		while(getline(openFile, ln)){
			if(!ln.length() || ln[0] == '#'){continue;}

			string fdir, fscf, sdir, sscf, rfscf, rsscf;
			vector<string> col = Split_by_space(ln);
			
			fscf = col[0];
			sscf = col[1];
			float score = atof(col[2].c_str());

			vector<string> reverse_crds = Return_reverse (fscf, sscf);
			rfscf = reverse_crds[0]; 
			rsscf = reverse_crds[1];

			Save_info(fscf, sscf, score);	
			Save_info(rfscf, rsscf, score);	
		}
		openFile.close();
	}
}

int _Network::Search_localgreedypath(const string node, const string goal, float score, map<vector<string>, float>& path){	
	if(node_list.find(node) != node_list.end()){
		node_list[node]->Update_visit(n_visit);
	}
	next_visit.push_back(node);
	
	if(goal == node){
		float cnt = next_visit.size();
		
		path[next_visit] = (score/cnt);
		if(! next_visit.empty()) next_visit.pop_back() ;

		return 1;
	}
	
	if(node_list.find(node) != node_list.end()){
		map<string, float> adjnodes = node_list[node]->Get_adjNode();
			
		vector<string> maxnodes; int maxscore = -1;
		for(map<string, float>::iterator it = adjnodes.begin(); it != adjnodes.end(); ++it){
			string curnode = it->first;
			float curscore = it->second;

			vector<string> rcurnode = Return_reverse ("NA", curnode);
			
			if(curscore > maxscore){
				maxscore = curscore;
				maxnodes.clear();
				maxnodes.push_back(curnode);
			}else if(curscore == maxscore){
				maxnodes.push_back(curnode);
			}
		}
		
		for(int i=0; i<maxnodes.size(); ++i){
			string adjnode, radjnode, cur_dir, cur_id;
			adjnode = maxnodes[i];
			
			vector<string> reverse_crds = Return_reverse(node, adjnode);
			radjnode = reverse_crds[0];
			score += maxscore;

			if(adjnode == goal || node_list.find(adjnode) == node_list.end() || (!node_list[adjnode]->Is_visit() && !node_list[radjnode]->Is_visit())){
				Search_localgreedypath (adjnode, goal, score, path);
				
				if(node_list.find(adjnode) != node_list.end()){
					node_list[adjnode]->Update_visit(n_unvisit);
				}
			}
			score -= maxscore;
		}
	}
	next_visit.pop_back();
}
	
map<string, vector<string> > _Network::Search_globalgreedypath (const vector<string>& bs, const vector<string>& es){
	vector<_Link> unused_links;
	map<string, string> used_links;
	map<string, vector<string> > paths;

	for(map<string, _Node*>::iterator mit=node_list.begin(); mit!=node_list.end(); ++mit){
		string nid = mit->first;
		_Node* node = mit->second;
		if(node->Is_used()){ continue; }
		
		map<string, float> adjnodes = node->Get_adjNode();
		for(map<string, float>::iterator mit=adjnodes.begin(); mit!=adjnodes.end(); ++mit){
			string adjnode = mit->first;
			float score = mit->second;

			vector<string> reverse_crds = Return_reverse(nid, adjnode);
			_Link check_link = { reverse_crds[0], reverse_crds[1], 0 };

			if((find(unused_links.begin(), unused_links.end(), check_link) != unused_links.end()) || ((find(bs.begin(), bs.end(), nid) != bs.end() || find(bs.begin(), bs.end(), reverse_crds[0]) != bs.end()) || (find(es.begin(), es.end(), adjnode) != es.end() || find(es.begin(), es.end(), reverse_crds[1]) != es.end()))){ 
				continue; 
			} 
			
			_Link new_link = { nid, adjnode, score };
			unused_links.push_back(new_link);
		}
	}

	sort(unused_links.begin(), unused_links.end());
	
	for(int i=0; i<unused_links.size(); ++i){
		vector<string> delete_ele;

		string src = unused_links[i].src, tar = unused_links[i].tar;
		vector<string> reverse_crds = Return_reverse (src, tar);
		string rsrc = reverse_crds[0], rtar = reverse_crds[1];

		if(used_links.find(src) != used_links.end() || used_links.find(rsrc) != used_links.end()){ continue; }
	
		string connect = "none";
		for(map<string, vector<string> >::iterator mit=paths.begin(); mit != paths.end(); ++mit){	
			string begin = mit->first;
			vector<string> cur_path = mit->second;
			string end = cur_path[cur_path.size()-1];
		
			if(end == src || end == rsrc){
				if(connect != "none"){
					if(end == src){
						if(tar == connect || src == connect){
							if(paths[begin][paths[begin].size()-1] == paths[connect][0]){
								paths[begin].insert(paths[begin].end(), paths[connect].begin()+1, paths[connect].end());		
							}else{
								paths[begin].insert(paths[begin].end(), paths[connect].begin(), paths[connect].end());		
							}
						}else{
							vector<string> reverse_path = paths[connect];
							reverse(reverse_path.begin(), reverse_path.end());
						
							for(int i=0; i<reverse_path.size(); ++i){
								vector<string> reverse_crds = Return_reverse ("NA", reverse_path[i]);
								reverse_path[i] = reverse_crds[0];
							}
	
							if(paths[begin][paths[begin].size()-1] == reverse_path[0]){
								 paths[begin].insert(paths[begin].end(), reverse_path.begin()+1, reverse_path.end());
							}else{
							 	paths[begin].insert(paths[begin].end(), reverse_path.begin(), reverse_path.end());
							}
						}
					}else{
						if(connect == rsrc || connect == rtar){
							if(paths[begin][paths[begin].size()-1] == paths[connect][0]){
								paths[begin].insert(paths[begin].end(), paths[connect].begin()+1, paths[connect].end());		
							}else{
								paths[begin].insert(paths[begin].end(), paths[connect].begin(), paths[connect].end());		
							}
						}else{
							vector<string> reverse_path = paths[connect];
							reverse(reverse_path.begin(), reverse_path.end());
						
							for(int i=0; i<reverse_path.size(); ++i){
								vector<string> reverse_crds = Return_reverse ("NA", reverse_path[i]);
								reverse_path[i] = reverse_crds[0];
							}
							
							if(paths[begin][paths[begin].size()-1] == reverse_path[0]){
								 paths[begin].insert(paths[begin].end(), reverse_path.begin()+1, reverse_path.end());
							}else{
							 	paths[begin].insert(paths[begin].end(), reverse_path.begin(), reverse_path.end());
							}
						}
					}

					delete_ele.push_back(connect);
				}else{
					if(end == src){
						paths[begin].push_back(tar);
					}else{
						paths[begin].push_back(rtar);
					}
				}
				
				connect = begin;
				used_links[src] = tar;
				used_links[rsrc] = rtar;
			}

			if(begin == tar || begin == rtar){
				if(connect != "none"){
					if(begin == tar){
						if(connect == rsrc || connect == rtar){
							vector<string> reverse_path = paths[connect];
							reverse(reverse_path.begin(), reverse_path.end());
						
							for(int i=0; i<reverse_path.size(); ++i){
								vector<string> reverse_crds = Return_reverse ("NA", reverse_path[i]);
								reverse_path[i] = reverse_crds[0];
							}
							
							delete_ele.push_back(connect);
							connect = reverse_path[0];
							paths[reverse_path[0]] = reverse_path;
								
							if(paths[reverse_path[0]][paths[reverse_path[0]].size()-1] == paths[begin][0]){
								paths[reverse_path[0]].insert(paths[reverse_path[0]].end(), paths[begin].begin()+1, paths[begin].end());
							}else{
								paths[reverse_path[0]].insert(paths[reverse_path[0]].end(), paths[begin].begin(), paths[begin].end());
							}
						}else{
							if(paths[connect][paths[connect].size()-1] == paths[begin][0]){
								paths[connect].insert(paths[connect].end(), paths[begin].begin()+1, paths[begin].end());
							}else{
								paths[connect].insert(paths[connect].end(), paths[begin].begin(), paths[begin].end());
							}
						}
					}else{
						if(connect == src || connect == tar){
							vector<string> reverse_path = paths[connect];
							reverse(reverse_path.begin(), reverse_path.end());
						
							for(int i=0; i<reverse_path.size(); ++i){
								vector<string> reverse_crds = Return_reverse ("NA", reverse_path[i]);
								reverse_path[i] = reverse_crds[0];
							}
							
							delete_ele.push_back(connect);
							connect = reverse_path[0];
							paths[reverse_path[0]] = reverse_path;

							if(reverse_path[reverse_path.size()-1] == begin){
								paths[reverse_path[0]].insert(paths[reverse_path[0]].end(), paths[begin].begin()+1, paths[begin].end());
							}else{
								paths[reverse_path[0]].insert(paths[reverse_path[0]].end(), paths[begin].begin(), paths[begin].end());
							}
						}else{
							if(paths[connect][paths[connect].size()-1] == begin){
								paths[connect].insert(paths[connect].end(), paths[begin].begin()+1, paths[begin].end());
							}else{
								paths[connect].insert(paths[connect].end(), paths[begin].begin(), paths[begin].end());
							}
						}
					}
					delete_ele.push_back(begin);
				}else{
					if(begin == tar){
						cur_path.insert(cur_path.begin(), src);
						paths[src] = cur_path;
						
						connect = src;
					}else{
						cur_path.insert(cur_path.begin(), rsrc);
						paths[rsrc] = cur_path;

						connect = rsrc;
					}
					delete_ele.push_back(begin);
				}
				
				used_links[src] = tar;
				used_links[rsrc] = rtar;
			}
		}

		if(delete_ele.size()){
			for(int i=0; i<delete_ele.size(); i++){
				paths.erase(delete_ele[i]);
			}
		}
		
		if(connect == "none"){
			vector<string> new_path { src, tar };
			paths[src] = new_path;
			
			used_links[src] = tar;
			used_links[rsrc] = rtar;
		}
	}
	
	map<string, vector<string>> reverse_path;
	for(map<string, vector<string> >::iterator mit=paths.begin(); mit != paths.end(); ++mit){	
		vector<string> path = mit->second;
		
		reverse(path.begin(), path.end());
		for(int i=0; i<path.size(); ++i){
			vector<string> reverse_crd = Return_reverse("NA", path[i]);
			path[i] = reverse_crd[0];
		}
		reverse_path[path[0]] = path;
	}
	
	paths.insert(reverse_path.begin(), reverse_path.end());

	return paths;
}

class _Traveler {
	map<string, string> anchor;
	vector<string> ends;
	vector<string> begins;

  public:
	int Get_anchor(string bid_f);
	int Fill_synteny_link (_Network& net);
	int Extend_synteny_link (_Network& net);
};

int _Traveler::Get_anchor(string bid_f){
	string ln, pre_chr = "", pre_id = "";
	map<string, string> id_convert;

	ifstream openFile1(bid_f.data());
	if(!openFile1.is_open()){
		cerr << "[Error] Cannot open " << bid_f << endl;
		return 0;
	}

	cerr << "Read file " << bid_f << endl;
	while(getline(openFile1, ln)){
		if(!ln.length() || ln[0] == '#'){ continue; }
		vector<string> col = Split_by_space(ln);
		
		for(int i=1; i<col.size()-1; ++i){
			string preid = col[i-1], id = col[i];
			anchor[preid] = id;
		}
		begins.push_back(col[0]);
		ends.push_back(col[col.size()-2]);
	}
	openFile1.close();
}

int _Traveler::Fill_synteny_link (_Network& net){
	string ln;

	for(map<string, string>::iterator mit = anchor.begin(); mit != anchor.end(); ++mit){
		string fscf = mit->first, sscf = mit->second;
		map<vector<string>, float> paths;
		
		string cur_dir, cur_id;
		
		vector<string> reverse_crds = Return_reverse(fscf, sscf);
		string rfscf = reverse_crds[0];

		if(net.Search_node(fscf) && net.Search_node(rfscf)){
			net.Search_localgreedypath(fscf, sscf, 0, paths);
		}
		
//		if(! paths.size()){ continue; }
		float max_score = -1;
		
		for(map<vector<string>, float>::iterator smit=paths.begin(); smit!=paths.end(); ++smit){
			vector<string> next_visit = smit->first;
			float score = smit->second;
			if(score > max_score) { max_score = score; }
		}

		vector<string> max_path;
		for(map<vector<string>, float>::iterator smit=paths.begin(); smit!=paths.end(); ++smit){
			vector<string> cur_path = smit->first;
			float score = smit->second;
			
			if(score == max_score && max_path.size() < cur_path.size()){ max_path = cur_path; }	
		}
	
		cout << "#" << fscf << "\t" << sscf << endl;
		for(int i=0; i<max_path.size(); ++i){
			if(net.node_list.find(max_path[i]) != net.node_list.end()){
				net.node_list[max_path[i]]->Update_used();
			}	
			
			vector<string> reverse_crds = Return_reverse ("NA", max_path[i]);
			if(net.node_list.find(reverse_crds[0]) != net.node_list.end()){
				net.node_list[reverse_crds[0]]->Update_used();
			}

			if(i){ cout << " "; } 
			cout <<  max_path[i];
		}
		cout << endl;
	}
}

int _Traveler::Extend_synteny_link (_Network& net){
	map<string, vector<string> > path_candidate = net.Search_globalgreedypath(begins, ends);

	for(int i=0; i<ends.size(); ++i){
		string cur_end = ends[i];

		vector<string> cur_path;
		if(path_candidate.find(cur_end) != path_candidate.end()){
			cur_path = path_candidate[cur_end];
			
			vector<string> reverse_crd = Return_reverse ("NA", cur_path[cur_path.size()-1]);
			string path_end = reverse_crd[0];

			path_candidate.erase(cur_end);
			path_candidate.erase(path_end);
		}

		if(cur_path.size() > 0){
			cout << "#" << ends[i] << " $" << endl;
			cout << cur_path[0];
			for(int j=1; j<cur_path.size()-1; ++j){
				cout << " " << cur_path[j];
			}
			cout << " $" << endl;
		}
	}

	for(int i=0; i<begins.size(); ++i){
		string cur_beg = begins[i];
		vector<string> reverse_crd = Return_reverse ("NA", cur_beg);
		string reverse_beg = reverse_crd[0];
	
		vector<string> cur_path;
		if(path_candidate.find(reverse_beg) != path_candidate.end()){
			cur_path = path_candidate[reverse_beg];
			
			vector<string> reverse_crd = Return_reverse ("NA", cur_path[cur_path.size()-1]);
			string path_end = reverse_crd[0];

			path_candidate.erase(reverse_beg);
			path_candidate.erase(path_end);
		}
		
		if(cur_path.size() > 0){
			cout << "#" << reverse_beg << " $" << endl;
			cout << cur_path[0];
			for(int j=1; j<cur_path.size(); ++j){
				cout << " " << cur_path[j];
			}
			cout << " $" << endl;
		}
	}
}

int main(int argc, char** argv){	
	if(argc != 6){
		cerr << "[Error] fill_extension perform_fill[1|0] perform_extension[1|0] bid_file SR_link LR_link" << endl;
		return 0;
	}
	
	string bid_file = argv[3], sr_link = argv[4], lr_link = argv[5];
	int perform_fill = atoi(argv[1]), perform_ext = atoi(argv[2]);

	vector<string> network_src;
	_Network scf_network;
	_Traveler walker;

	if(sr_link != "NA"){
		network_src.push_back(sr_link);
	}

	if(lr_link != "NA"){ 
		network_src.push_back(lr_link);
	}

	scf_network.Build_network(network_src);
//	scf_network.Print_network();

	walker.Get_anchor(bid_file);
	if(perform_fill){
		walker.Fill_synteny_link (scf_network);
	}
	if(perform_ext){
		walker.Extend_synteny_link (scf_network);
	}
}

vector<string> Split_by_space (const string& str){
	vector<string> tokens;
	stringstream sst;
	sst.str(str);
	
	string token;
	while(sst >> token){
		tokens.push_back(token);
	}
	return tokens;
}

vector<string> Return_reverse (string crd_a, string crd_b){
	vector<string> reverse;

	smatch match;
	regex re("(-*)(\\S+)");
	string fdir, fscf, sdir, sscf, rfscf, rsscf;

	if(regex_match(crd_a, match, re)){
		fdir = match[1];
		fscf = match[2];
	}
	if(regex_match(crd_b, match, re)){
		sdir = match[1];
		sscf = match[2];
	}
	
	rfscf = (sdir == "-")? sscf : "-"+sscf;
	rsscf = (fdir == "-")? fscf : "-"+fscf;

	reverse.push_back(rfscf);
	reverse.push_back(rsscf);

	return reverse;
}
