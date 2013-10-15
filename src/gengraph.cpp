#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <string.h> 
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <stdarg.h>
#include <time.h>
#include <vector>
#include <gvc.h>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>
#include <boost/graph/connected_components.hpp>

#include "springlayout.hpp"
#include "smithwaterman.hpp"

using namespace alignment;
using namespace boost;
using namespace std;

#define VERSION 0.1
#define SEQ_ID_THRESHOLD 90.0

typedef boost::square_topology<>::point_type Point;
struct VertexProperties{
    size_t index, topology, hom_group;
    string name, location;
    double x, y;
    Point point;
};

typedef adjacency_list< vecS, vecS, undirectedS, VertexProperties, property< edge_weight_t, double > > Graph;
typedef boost::property_map<Graph, std::size_t VertexProperties::*>::type VertexIndexPropertyMap;
typedef boost::property_map<Graph, std::string VertexProperties::*>::type VertexNamePropertyMap;
typedef boost::property_map<Graph, Point VertexProperties::*>::type PositionMap;
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightPropertyMap;
typedef boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;

struct Residue{
	Residue(string r, int n, double x, double y, double z) : res_(r), num_(n), x_(x), y_(y), z_(z) {}
	string res_;
	int num_;
	double x_, y_, z_;
};

class PDB{

public:
	PDB();
	explicit PDB(string, double s = 35.0);
	~PDB();
	void init();
	void parse_pdb(string, double s = 35.0, int tm_only = 0);
	void build_graph(const double dist = 8.0, const int threshold = 10, const double scale = 0.075);
	int layout_graph(GVC_t*, bool);
	void write_graph(GVC_t*);
	void write_subgraph(GVC_t*,const string&);
	void find_symmetry(int search_type = 1, int sym_max = 27, int tm_only = 0);
	void print_graph(const double radius);
	void free_layout(GVC_t*);
	void print_chains();
	void circle_layout();
	void kamada_layout();
	void set_tm_chains(const string&, int);
	int get_chain_count() const;
	vector<string> get_chains() const;
	map<string,int> tm_chains;
	map<string,int> homology_cluster_id;
	vector<int> correspond;
	string base, pdbname, output;
	Graph g;
	Agraph_t* gv;	
	vector<Agnode_t*> gv_nodes;
	VertexIndexPropertyMap vertexHomologyGroupMap;
	int layout;

private:
	VertexIndexPropertyMap vertexIdPropertyMap;
	VertexIndexPropertyMap vertexTopologyPropertyMap;
	VertexNamePropertyMap vertexNamePropertyMap;
	VertexNamePropertyMap vertexLocationPropertyMap;	
	vector<string> chains;
	map <string, string> aa3_1, aa1_3; 
	map <int, string> num_to_chain;
	map<string,string> chain_to_seq;	
	map<string,int> chain_lengths, chain_to_num;	
	map<string,vector<Residue> > chain_to_ca;
	dynamic_properties dp;		
	Agnode_t *n;
	Agedge_t *e;
	int hom_clusters;
};

PDB::PDB(){

	init();
}

PDB::PDB(string target, double s){

	init();
	parse_pdb(target,s);
}

PDB::~PDB(){
	//cout << "Destroying stuff..." << endl;
	agclose(gv);
}

void PDB::free_layout(GVC_t* gvc){

	// Despite this call we still leak some memory
	// apparently caused by some 3rd party libraries
	// http://www.graphviz.org/mantisbt/view.php?id=2356
	gvFreeLayout(gvc, gv);  

}

vector<string> PDB::get_chains() const {

	return(chains);
}

void PDB::set_tm_chains(const string& t, int c){

	tm_chains[t] = c;
}

int PDB::get_chain_count() const {

	return(chains.size());
}

void PDB::init(){

	static const char* aa3[] = {"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};
	static const char* aa1[] = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};
	for(int d = 0; d < 20; d++){
		aa3_1[aa3[d]] = aa1[d];
		aa1_3[aa1[d]] = aa3[d];
	}
	// Graphviz	
	gv = agopen((char*)base.c_str(), Agundirected, 0);
	agsafeset(gv, (char*)"outputorder", (char*)"edgesfirst", (char*)"");
	agsafeset(gv, (char*)"splines", (char*)"true", (char*)"");
	//agsafeset(gv, (char*)"size", (char*)"20,20", (char*)"");
	//agsafeset(gv, (char*)"ratio", (char*)"fill", (char*)"");
	//agsafeset(gv, (char*)"overlap", (char*)"false", (char*)"");
	//agsafeset(gv, (char*)"overlap", (char*)"scale", (char*)"");
	//agsafeset(gv, (char*)"overlap", (char*)"scalexy", (char*)"");
	//agsafeset(gv, (char*)"overlap", (char*)"orthoxy", (char*)"");
	layout = 1;
	vertexIdPropertyMap = boost::get(&VertexProperties::index, g);
	vertexTopologyPropertyMap = boost::get(&VertexProperties::topology, g);
	vertexNamePropertyMap = boost::get(&VertexProperties::name, g);
	vertexLocationPropertyMap = boost::get(&VertexProperties::location, g);	
	vertexHomologyGroupMap = boost::get(&VertexProperties::hom_group, g);
}

void PDB::parse_pdb(string target, double seq_thresh, int tm_only){

	cout << "Parsing " << target << ":" << endl;
	base = target;
	string empty = "";	
	base = regex_replace(base, boost::regex("\\.pdb$",boost::regex::perl|boost::regex::icase), empty);
	map<string,int>::const_iterator it;
	boost::filesystem::path pdbpath(base);
	pdbname = pdbpath.filename().string();
	char buf[1024];

	// Parse PDB file
	FILE *finpdb = fopen(target.c_str(), "r");
	if (finpdb != NULL){
		while (fgets (buf, 1024, finpdb)){
			string line = buf;
			if (strncmp(line.substr(0,4).c_str(),"ATOM",4) == 0){

				string atom = line.substr(13,4);
				erase_all(atom, " ");
				string res = line.substr(17,3);
				string chain = line.substr(21,1);

				it = tm_chains.find(chain);
				if(it != tm_chains.end() || (it == tm_chains.end() && !tm_only)){

					chains.push_back(chain);
					if(aa3_1.find(res) == aa3_1.end()){
						cout << "Unknown residue  " << res << ":" << endl << line;	
						res = "X";
					}

					if((strcmp (atom.c_str(),"CA") == 0 && strcmp (res.c_str(),"GLY") == 0)||(strcmp (atom.c_str(),"CB") == 0 && strcmp (res.c_str(),"GLY") != 0)){
						double x = atof(line.substr(30,8).c_str());
						double y = atof(line.substr(38,8).c_str());
						double z = atof(line.substr(46,8).c_str());	
						string resnum = line.substr(22,4);
						erase_all(resnum, " ");
						int n = boost::lexical_cast<int>(resnum);
						Residue r(aa3_1[res],n,x,y,z);

						if(chain_to_seq.find(chain) == chain_to_seq.end()){
							chain_to_seq[chain]	= aa3_1[res];
							chain_lengths[chain] = 1;							
							vector<Residue> tmp;
							tmp.push_back(r);
							chain_to_ca.insert(pair<string,vector<Residue> >(chain,tmp));
						}else{
							chain_to_seq[chain].append(aa3_1[res]);
							chain_lengths[chain]++;
							chain_to_ca[chain].push_back(r);	
						}	
					}

				}				
			}			
		}	
		fclose(finpdb);		
	}else{
		cout << "Couldn't open file " << target << endl;
	}	
	chains.erase(unique(chains.begin(),chains.end()),chains.end());

	// Detect homologous chains and assign them to a cluster
	// Each cluster will then have its own node shape
	int shape = 0;
	map<string,double> homology_cluster_pct;
	map<string,double>::const_iterator itd;
	for(unsigned int i = 0; i < chains.size(); i++){
		for(unsigned int j = i+1; j < chains.size(); j++){

			// Make sure localisation is the same too
			// We could also make sure the number of TM helices is the same...
			bool loc1 = false, loc2 = false;
			it = tm_chains.find(chains[i]);
			if(it != tm_chains.end()){
				loc1 = true;
			}
			it = tm_chains.find(chains[j]);
			if(it != tm_chains.end()){
				loc2 = true;
			}

			// Do Smith-Waterman alignemnt
			double pct_id = smithwaterman(chain_to_seq[chains[i]],chain_to_seq[chains[j]],2.0,-0.5,-1.0,false);
			if(100*pct_id > seq_thresh && loc1 == loc2){
				//cout << "Chains " << chains[i] << " and " << chains[j] << " are homolgous:\t\t";
				//printf("%3.2f %% ID\n",100*pct_id);				

				itd = homology_cluster_pct.find(chains[j]);
				if(itd == homology_cluster_pct.end()){
					homology_cluster_pct[chains[j]] = pct_id;
				}	

				it = homology_cluster_id.find(chains[i]);
				if(it == homology_cluster_id.end()){					
					//cout << chains[i] << " is a new chain, shape is now " << shape << endl;
					//cout << chain_to_seq[chains[i]] << endl;
					homology_cluster_id[chains[j]] = homology_cluster_id[chains[i]] = shape;
					//cout << "Adding chain " << chains[j] << " to homology map with ID " << pct_id << endl;
					//homology_cluster_pct[chains[j]] = pct_id;
					shape++;
					hom_clusters++;
					// Out of shapes
					if(shape == 17) shape = 0;
				}else{
					if(pct_id >= homology_cluster_pct[chains[j]]){
						//cout << "Seq. ID is >= " << pct_id << " vs " << homology_cluster_pct[chains[j]] << " - changing homology group." << endl;	
						homology_cluster_pct[chains[j]] = pct_id;
						//cout << "Setting chain " << chains[j] << " to same node type as chain " << chains[i] << " (" << homology_cluster_id[chains[i]] << ")" << endl;
						homology_cluster_id[chains[j]] = homology_cluster_id[chains[i]];							
					}else{
						//cout << "Seq. ID is lower " << pct_id << " vs " << homology_cluster_pct[chains[j]] << " - NOT changing homology group." << endl;	
					}

				}
			}else{
				//cout << "Chains " << chains[i] << " and " << chains[j] << " are not homolgous:\t";
				//printf("%3.2f %% ID\n",100*pct_id);
				//if(loc1 != loc2) cout << "(or opposite localisation)" << endl;
			}
			//printf("%3.2f %% ID\n",100*pct_id);
		}
		it = homology_cluster_id.find(chains[i]);
		if(it == homology_cluster_id.end()){
			homology_cluster_id[chains[i]] = shape;
			shape++;
			hom_clusters++;
			// Out of shapes
			if(shape == 17) shape = 0;
		}
	}
	
	cout << endl;
	for(unsigned int i = 0; i < chains.size(); i++){
		cout << "Chain " << chains[i] << " ---> " << homology_cluster_id[chains[i]] << endl;
	}
	cout << hom_clusters << " homologous clusters found." << endl;
	cout << endl;
}

double get_dist(Residue& a, Residue& b){
	double dx = a.x_ - b.x_;
	double dy = a.y_ - b.y_;
	double dz = a.z_ - b.z_;
	return(sqrt(dx*dx + dy*dy + dz*dz));

}

void PDB::build_graph(const double dist, const int threshold, const double scale){

	map<string,int>::const_iterator it;

	static const char* shapes[] = {"circle","square","house","pentagon","hexagon","octagon","septagon","diamond","ellipse","oval","egg","triangle","invtriangle","invtrapezium","invhouse","trapezium","parallelogram"};

	for(int i = 0; i < get_chain_count(); i++){
		chain_to_num[get_chains()[i]] = i;
		num_to_chain[i] = get_chains()[i];
        VertexDescriptor vd = boost::add_vertex(g);
        vertexIdPropertyMap[vd] = i;
        vertexNamePropertyMap[vd] = get_chains()[i];

        // Graphviz
        n = agnode(gv, (char*)get_chains()[i].c_str(), 1);
        agsafeset(n, (char*)"name", (char*)get_chains()[i].c_str(), (char*)"");

        //cout << shapes[homology_cluster_id[get_chains()[i]]] << endl;

        agsafeset(n, (char*)"shape", (char*)shapes[homology_cluster_id[get_chains()[i]]], (char*)"");
        //agsafeset(n, (char*)"shape", (char*)"circle", (char*)"");
        agsafeset(n, (char*)"fontcolor", (char*)"white", (char*)"");
        agsafeset(n, (char*)"fontname", (char*)"helvetica", (char*)"");
        agsafeset(n, (char*)"color", (char*)"grey", (char*)"");
        agsafeset(n, (char*)"fixedsize", (char*)"true", (char*)"");

        vertexHomologyGroupMap[vd] = homology_cluster_id[get_chains()[i]];

        if(tm_chains.find(get_chains()[i]) != tm_chains.end()){
        	vertexLocationPropertyMap[vd] = "membrane";
        	vertexTopologyPropertyMap[vd] = tm_chains[get_chains()[i]];
        	agsafeset(n, (char*)"style", (char*)"radial", (char*)"");
        	agsafeset(n, (char*)"fillcolor", (char*)"mediumblue;0.1:steelblue;0.9", (char*)"");
        	string size = to_string(0.5+scale*tm_chains[get_chains()[i]]);
        	agsafeset(n, (char*)"width", (char*)size.c_str(), (char*)"");
        	agsafeset(n, (char*)"height", (char*)size.c_str(), (char*)"");
        }else{
        	vertexLocationPropertyMap[vd] = "globular";	
        	vertexTopologyPropertyMap[vd] = 0;
        	agsafeset(n, (char*)"style", (char*)"radial", (char*)"");
        	agsafeset(n, (char*)"fillcolor", (char*)"red;0.9:firebrick;0.1", (char*)"");
        	agsafeset(n, (char*)"width", (char*)to_string(0.5).c_str(), (char*)"");
        	agsafeset(n, (char*)"height", (char*)to_string(0.5).c_str(), (char*)"");
        }
        gv_nodes.push_back(n);
	}
	cout << num_vertices(g) << " vertices added." << endl;
	
    for (map<string,int>::iterator it = chain_lengths.begin() ; it != chain_lengths.end(); ++it){
    	for (map<string,int>::iterator it2 = it ; it2 != chain_lengths.end(); ++it2){
    		if(it != it2){
    			vector<int> contacts_c1, contacts_c2;    			
		    	BOOST_FOREACH(Residue & r, chain_to_ca[it->first]){
		    		BOOST_FOREACH(Residue & q, chain_to_ca[it2->first]){
		    			if(get_dist(r,q) < dist){
		    				contacts_c1.push_back(r.num_);	
		    				contacts_c2.push_back(q.num_);	
		    			}	
		    		}
		    	}		    	
		    	contacts_c1.erase(unique(contacts_c1.begin(),contacts_c1.end() ),contacts_c1.end());
		    	contacts_c2.erase(unique(contacts_c2.begin(),contacts_c2.end() ),contacts_c2.end());
		    	int c = contacts_c1.size() + contacts_c2.size();
		    	if(c >= threshold){
		    		//cout << chain_to_num[it->first] << " vs " << chain_to_num[it2->first] << " contact (" << c << ")" << endl;
		    		add_edge(vertex(chain_to_num[it->first],g), vertex(chain_to_num[it2->first],g), 1.0, g);

		    		// Graphviz
		    		e = agedge(gv,gv_nodes[chain_to_num[it->first]],gv_nodes[chain_to_num[it2->first]],0,1);
		    		agsafeset(e, (char*)"color", (char*)"grey", (char*)"");
		    		// weight not used by sfpd
		    		agsafeset(n, (char*)"weight", (char*)to_string(2.8).c_str(), (char*)"");
		    	}
    		}
    	}	
    }    
    cout << num_edges(g) << " edges added." << endl;

    //if(!num_vertices(g)) return;
}

void PDB::circle_layout(){

	//cout << "In circle..." << endl;
	PositionMap positionMap = boost::get(&VertexProperties::point, g);

    boost::circle_graph_layout(g, positionMap, 300);
    
	int x = 0;
	boost::graph_traits<Graph>::vertex_iterator i, end;
	for(boost::tie(i, end) = boost::vertices(g); i != end; ++i) {
		//cout << positionMap[*i][0] << "," << positionMap[*i][1] << endl;
		double n = positionMap[*i][0];
		// Alter denominator here to prevent node overlaps
		//n = static_cast<double>(static_cast<int>(n*1000))/100;
		double m = positionMap[*i][1];
		//n += 400;
		//m = static_cast<double>(static_cast<int>(m*1000))/100;

		string pos = to_string(n) + "," + to_string(m);
		//cout << pos << endl;
		agsafeset(gv_nodes[x], (char*)"pos", (char*)pos.c_str(), (char*)"");
		x++;
		put(&VertexProperties::x, g, *i, positionMap[*i][0]*10);
		put(&VertexProperties::y, g, *i, positionMap[*i][1]*10);
	}
	//cout << "done..." << endl;
	
}	

void PDB::kamada_layout(){

    PositionMap positionMap = boost::get(&VertexProperties::point, g);
    WeightPropertyMap weightPropertyMap = boost::get(edge_weight_t(), g);

    boost::circle_graph_layout(g, positionMap, 100);
    double tolerance = 0.00001;
    int iterations = 100000;
    layout_and_iteration_tolerance<double> tol(tolerance,iterations);
    //cout << "Calling Kamada Kawai spring layout..." << endl;
    bool retval = boost::kamada_kawai_spring_layout(g, positionMap, weightPropertyMap, 
                                                       boost::square_topology<>(), boost::side_length<double>(10),     
                                                       tol, 1, vertexIdPropertyMap);
    if (!retval){
        cout << "kamada_kawai_spring_layout returned false!" << endl;
        return;
    }
    
	int x = 0;
	boost::graph_traits<Graph>::vertex_iterator i, end;
	for(boost::tie(i, end) = boost::vertices(g); i != end; ++i) {
		//cout << positionMap[*i][0] << "," << positionMap[*i][1] << endl;
		double n = positionMap[*i][0];
		// Alter denominator here to prevent node overlaps
		n = static_cast<double>(static_cast<int>(n*1000))/20;
		double m = positionMap[*i][1];
		n += 400;
		m = static_cast<double>(static_cast<int>(m*1000))/20;
		string pos = to_string(n) + "," + to_string(m);
		agsafeset(gv_nodes[x], (char*)"pos", (char*)pos.c_str(), (char*)"");
		x++;
		put(&VertexProperties::x, g, *i, positionMap[*i][0]*10);
		put(&VertexProperties::y, g, *i, positionMap[*i][1]*10);
	}
	//cout << "Finished KK." << endl;
	// http://levelfour.googlecode.com/svn/branches/dev/vd2/Components/D2M/sandratest/source/GraphPropertyMapper.h
  	// http://levelfour.googlecode.com/svn/branches/dev/vd2/Components/D2M/sandratest/source/TestPlotter.h
    // https://code.google.com/p/levelfour/source/browse/branches/dev/vd2/Components/D2M/sandratest/kamada_layout_example.cpp
  
    dp.property("Label", boost::get(&VertexProperties::name, g));
    dp.property("Location", boost::get(&VertexProperties::location, g));
    dp.property("Topology", boost::get(&VertexProperties::topology, g));
    dp.property("X", boost::get(&VertexProperties::x, g));
    dp.property("Y", boost::get(&VertexProperties::y, g));
    dp.property("Weight", boost::get(edge_weight_t(), g));   
	
}	

int PDB::layout_graph(GVC_t* gvc, bool t){

    if(!num_edges(g)){
    	cout << pdbname << " has no edges - skipping layout." << endl;
    	return(0);
    }
    if(num_vertices(g) < 2){
    	cout << pdbname << " has less than two vertices - skipping layout." << endl;
    	return(0);
    }
	
	//circle_layout();
    //gvLayout(gvc, gv, "nop");
    //gvFreeLayout(gvc, gv);  
    //return(1);

    switch(layout){
    	case 2:
    		if (t) cout << "Calling fdp layout..." << endl;
			//circle_layout();
		    //gvLayout(gvc, gv, "nop");
		    //gvFreeLayout(gvc, gv);  
    		gvLayout (gvc, gv, "fdp");
    		break;
    	case 3:
    		if (t) cout << "Calling circo layout..." << endl;
    		gvLayout (gvc, gv, "circo");
    		break;    	
    	case 4:
    		if (t) cout << "Calling neato layout..." << endl;

    		// Better results if arranged on circle first
			circle_layout();
		    gvLayout(gvc, gv, "nop");
		    gvFreeLayout(gvc, gv);  
    		gvLayout (gvc, gv, "neato");
    		break;   
    	case 5:
    		if (t) cout << "Calling kamada-kawai layout..." << endl;
    		kamada_layout();
    	 	gvLayout(gvc, gv, "nop");
    		break;   
    	default:
    		if (t) cout << "Calling sfdp layout..." << endl;
    		gvLayout (gvc, gv, "sfdp");
    }
    //cout << "Layout done." << endl;
    return(1);
}

void PDB::write_graph(GVC_t* gvc){

	// Boost/Graphml
	/*
    string out = base + ".graphml";
    cout << "Writing " << out << endl;  
    filebuf fb;
    fb.open (out.c_str(),ios::out);
    std::ostream os(&fb);
    write_graphml(os, g, dp, true);
    fb.close();
	*/

	string gv_out = base + ".svg";
	if(!output.empty()){
		if(output.substr(output.length()-1,1) != "/"){
			output += "/";
		}
		gv_out = output + pdbname + ".svg";
	}
	cout << "Writing " << gv_out << endl;
	gvRenderFilename (gvc, gv, "svg", (char*)gv_out.c_str());

	//gv_out = base + ".dot";
	//cout << "Writing " << gv_out << endl;
	//gvRenderFilename (gvc, gv, "dot", (char*)gv_out.c_str());

	/*
    Agnode_t *v;
    //gvRender (gvc, gv, "dot", NULL);
	for(v = agfstnode(gv); v; v = agnxtnode(gv,v)){
		char* tmp1 = agget(v,(char*)"pos");
		printf("X value is %s\n",tmp1);
	}
	*/	
	/*
	gv_out = base + ".png";
	cout << "Writing " << gv_out << endl << endl;
	gvRenderFilename (gvc, gv, "png", (char*)gv_out.c_str());
	*/
}

void PDB::write_subgraph(GVC_t* gvc, const string& pdb2){

	//cout << "MCS nodes: " << agnnodes(gv) << endl;
	//cout << "MCS gv_nodes[]: " << gv_nodes.size() << endl;
	//cout << "MCS edges: " << agnedges(gv) << endl << endl;
    // Delete nodes and edges not present in the MCS
	vector<int>::const_iterator it;
	for(unsigned int i = 0; i < gv_nodes.size(); i++){
		it = find(correspond.begin(),correspond.end(),i);
		if(it == correspond.end()){
			//cout << "deleting node " << i << endl;
			agdelnode(gv,gv_nodes[i]);
		}else{
			//cout << "not deleting node " << i << endl;
		}
	}
	//cout << "MCS nodes: " << agnnodes(gv) << endl << endl;
	//cout << "MCS gv_nodes[]: " << gv_nodes.size() << endl;
	//cout << "MCS edges: " << agnedges(gv) << endl;

    if(!agnedges(gv)){
    	cout << "MCS has no edges - skipping layout." << endl;
    	return;
    }
    if(agnnodes(gv) < 2){
    	cout << "MCS  has less than two vertices - skipping layout." << endl;
    	return;
    }
	string gv_out = base + "-" + pdb2 + "_mcs.svg";
	if(!output.empty()){
		if(output.substr(output.length()-1,1) != "/"){
			output += "/";
		}
		gv_out = output + pdbname + "-" + pdb2 + "_mcs.svg";
	}
	cout << "Writing " << gv_out << endl;
	gvRenderFilename (gvc, gv, "svg", (char*)gv_out.c_str());
	//gv_out = output + pdbname + "-" + pdb2 + "_mcs.dot";
	//gvRenderFilename (gvc, gv, "dot", (char*)gv_out.c_str());
}

void PDB::print_chains(){

    for (map<string,int>::iterator it = chain_lengths.begin() ; it != chain_lengths.end(); ++it){
    	cout << "Chain:\t" << it->first << endl;
    	cout << "Length:\t" << it->second << endl;
    	cout << "Seq:\t" << chain_to_seq[it->first] << endl << endl;
    	BOOST_FOREACH(Residue & r, chain_to_ca[it->first]){
    		cout << r.x_ << "\t" << r.y_ << "\t" << r.z_ << endl;
    	}
    	cout << endl;	
    }	
}

template <typename Iterator> inline bool next_combination(const Iterator first, Iterator k, const Iterator last){
	
	/* Credits: Thomas Draper */
	// http://stackoverflow.com/questions/5095407/n-choose-k-implementation
	if ((first == last) || (first == k) || (last == k))
		return false;

	Iterator itr1 = first;
	Iterator itr2 = last;
	++itr1;
	if (last == itr1)
		return false;
	itr1 = last;
	--itr1;
	itr1 = k;
	--itr2;
	while (first != itr1){
		if (*--itr1 < *itr2){
			Iterator j = k;
			while (!(*itr1 < *j)) ++j;
			iter_swap(itr1,j);
			++itr1;
			++j;
			itr2 = k;
			rotate(itr1,j,last);
			while (last != j){
				++j;
				++itr2;
			}
			rotate(k,itr2,last);
			return true;
		}
	}
	rotate(first,k,last);
	return false;
}

template <typename GraphFirst,typename GraphSecond> struct sym_callback {

	sym_callback(const GraphFirst& graph1,const GraphSecond& graph2, int* p, int* m, int* c1, int* c2) :
    	m_graph1(graph1), m_graph2(graph2){ 
		q = p;
		n = m;
		correspond1 = c1;
		correspond2 = c2;
    }

	template <typename CorrespondenceMapFirstToSecond, typename CorrespondenceMapSecondToFirst>
  			bool operator()(CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
            CorrespondenceMapSecondToFirst correspondence_map_2_to_1,
            typename graph_traits<GraphFirst>::vertices_size_type subgraph_size){

		// Print out correspondences between vertices	
		*q = subgraph_size;
		int ok = 1;
		BGL_FORALL_VERTICES_T(vertex1, m_graph1, GraphFirst) {
		  	if (get(correspondence_map_1_to_2, vertex1) != graph_traits<GraphSecond>::null_vertex()){
		   		//cout << vertex1 << " (" << correspond1[vertex1] << ") ----> " << get(correspondence_map_1_to_2, vertex1) << " (" << correspond2[get(correspondence_map_1_to_2, vertex1)] << ")" <<  endl;
				if(correspond1[vertex1] != correspond2[get(correspondence_map_1_to_2, vertex1)]){
					ok = 0;
					//cout << "FAIL!" << endl;
				}
			}	
		}
		if(ok){
		  	(*n)++;
		}		
		return (true);
  	}

private:
	int* q;
	int* n;
	int* correspond1;
	int* correspond2;	
    const GraphFirst& m_graph1;
    const GraphSecond& m_graph2;
};

struct simple_sort {
	bool operator() (int i,int j) {return (i<j);}
};

void PDB::find_symmetry(int search_type, int sym_max, int t){

	if(num_vertices(g) == 1){
		cout << pdbname << " only contains one chain - skipping symmetry search." << endl << endl;
		return;
	}
	if((int)num_vertices(g) > sym_max){
		cout << pdbname << " contains " << num_vertices(g) << " chains - skipping symmetry search." << endl << endl;
		return;
	}

	if(t){
		cout << "Searching for symmetry in " << pdbname << "..." << endl;
	}else{
		cout << "Searching for symmetry..." << endl;
	}

	int max_fold = 0;	
	
	int* c1 = new int[num_vertices(g)];
	for(unsigned int i = 0; i < num_vertices(g); i++){
		c1[i] = homology_cluster_id[get_chains()[i]];
	}

	graph_traits<Graph>::edge_iterator ei, ei_end;
	typedef graph_traits<Graph>::vertex_iterator vertex_iter;
	std::pair<vertex_iter, vertex_iter> vp;

	vector<int> ints(num_vertices(g));	
	//std::iota(ints.begin(), ints.end(), 0);
	for(unsigned int i = 0; i < num_vertices(g);i++){
		ints[i] = i;
	}

	//cout << "Homologous cluster size:" << endl;
	map<int,int> hom_clst_size;
	for(int i = 0; i < hom_clusters; i++){
		hom_clst_size[i] = 0;
	}
	for(unsigned int i = 0; i < num_vertices(g); i++){
		hom_clst_size[homology_cluster_id[get_chains()[i]]]++;
	}

	for(unsigned int k = num_vertices(g)/2; k >= 1; k--){
		
		if(num_vertices(g) % k == 0){
			//cout << "trying k = " << k << endl;
			vector<vector<int> > subunit_combinations;
			do{

				Graph g_s;
				VertexIndexPropertyMap vertexHomologyGroupMap_gs;
				vertexHomologyGroupMap_gs = boost::get(&VertexProperties::hom_group, g_s);
				vector<int> tmp;
				map<int,int> g_to_g_s,g_s_to_g;
				map<int,int> hom_clst_size_g_s;
				for(int i = 0; i < hom_clusters; i++){
					hom_clst_size_g_s[i] = 0;
				}

				// Add vertices
				for (unsigned int i = 0; i < k; ++i){
					//cout << get_chains()[ints[i]] << " ";
					tmp.push_back(ints[i]);
					//tmp2.push_back(homology_cluster_id[get_chains()[ints[i]]]);
					hom_clst_size_g_s[homology_cluster_id[get_chains()[ints[i]]]]++;
					//cout << homology_cluster_id[get_chains()[ints[i]]] << endl;
			        VertexDescriptor vd = boost::add_vertex(g_s);
					vertexHomologyGroupMap_gs[vd] = homology_cluster_id[get_chains()[ints[i]]];

			        g_to_g_s[ints[i]] = i;
			        g_s_to_g[i] = ints[i];
				}
				//cout << endl;
				
				// Check cluster sizes of subgraph are factors of cluster sizes of main graph
				int ok = 1;
				for(int i = 0; i < hom_clusters; i++){
					if(hom_clst_size_g_s[i]*(int)(num_vertices(g)/k) != hom_clst_size[i]){
						ok = 0;	
					}
				}
				
				if(ok){

					// Add edges between these nodes only
					for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
						if(find(tmp.begin(),tmp.end(), source(*ei, g)) != tmp.end() && find(tmp.begin(),tmp.end(), target(*ei, g)) != tmp.end()){
	  						add_edge(vertex(g_to_g_s[source(*ei, g)],g_s), vertex(g_to_g_s[target(*ei, g)],g_s), 1.0, g_s);
	  					}
	 				}	

				    vector<int> component(num_vertices(g_s));
	 				int num = connected_components(g_s, &component[0]);	 				
	 				
	 				if(num == 1 && search_type == 1){
	 					//cout << "1 component" << endl;
		 				int* p = new int(0);
						int* m = new int(0);
						int* s = new int(10000);						
						int* c2 = new int[num_vertices(g_s)];

						for(unsigned int i = 0; i < num_vertices(g_s); i++){
							c2[i] = homology_cluster_id[get_chains()[g_s_to_g[i]]];
						}

						sym_callback<Graph, Graph> my_callback(g, g_s, p, m, c1, c2);
						mcgregor_common_subgraphs_maximum_unique(g, g_s, true, my_callback,vertices_equivalent(make_property_map_equivalent(vertexHomologyGroupMap,vertexHomologyGroupMap_gs)));
						//mcgregor_common_subgraphs_maximum_unique(g, g_s, true, my_callback);

						if((*p == (int)k && ((unsigned)*m >= num_vertices(g)/k)) || num_vertices(g_s) == 1){
							//cout << num_vertices(g)/k << " fold symmetry: ";
							//cout << *m << " common subgraphs found with " << *p << " vertices" << endl;
							vector<int> tmp2;
							for (unsigned int i = 0; i < k; ++i){
								//cout << num_to_chain[ints[i]] << " ";
								tmp2.push_back(ints[i]);
							}	
							subunit_combinations.push_back(tmp2);
						}
						delete m;
						delete p;
						delete s;
						delete [] c2;
					}
					

					if(num == 1 && search_type != 1){
						vector<int> tmp2;
						for(unsigned int i = 0; i < k; ++i){
							//cout << num_to_chain[ints[i]] << " ";
							tmp2.push_back(ints[i]);
						}	
						subunit_combinations.push_back(tmp2);							
					}

				}

			}while(next_combination(ints.begin(),ints.begin() + k,ints.end()));

				// Combine subunits
				if(subunit_combinations.size()){
					do{
						vector<int> tmp;
						for (unsigned int i = 0; i < (num_vertices(g)/k); ++i){
							for(unsigned int j = 0; j < subunit_combinations[i].size(); j++){
								tmp.push_back(subunit_combinations[i][j]);
							}
						}
						sort(tmp.begin(), tmp.end(), simple_sort());
						tmp.erase(unique(tmp.begin(),tmp.end() ),tmp.end());
						if(tmp.size() == num_vertices(g)){
							cout << num_vertices(g)/k << " fold symmetry:" << endl;
							if((int)(num_vertices(g)/k) > max_fold){
								max_fold = 	num_vertices(g)/k;
							}
							for (unsigned int i = 0; i < (num_vertices(g)/k); ++i){
								cout << i+1 << " : ";
								for(unsigned int j = 0; j < subunit_combinations[i].size(); j++){
									cout << num_to_chain[subunit_combinations[i][j]] << " ";
								}
								cout << endl;
							}
						}
					}while(next_combination(subunit_combinations.begin(),subunit_combinations.begin() + (num_vertices(g)/k),subunit_combinations.end()));
				}
		}
	}
	if(!max_fold){
		cout << "No symmetry detected." << endl;
	}
	cout << endl;
	delete [] c1;
}

template <typename GraphFirst,typename GraphSecond> struct print_callback {

	print_callback(PDB* p1, PDB* p2, int* p, int* m, int* c1, int* c2, int* s) :
		prot1(p1), prot2(p2), m_graph1(prot1->g), m_graph2(prot2->g){
			q = p;
			n = m;
			correspond1 = c1;
			correspond2 = c2;
			init_score = s;
	}

	template <typename CorrespondenceMapFirstToSecond, typename CorrespondenceMapSecondToFirst>
  			bool operator()(CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
            CorrespondenceMapSecondToFirst correspondence_map_2_to_1,
            typename graph_traits<GraphFirst>::vertices_size_type subgraph_size){

		// Print out correspondences between vertices	
  		(*n)++;
  		int t = 0;
  		int score = 0;
    	BGL_FORALL_VERTICES_T(vertex1, m_graph1, GraphFirst) {
			// Skip unmapped vertices
    		if (get(correspondence_map_1_to_2, vertex1) != graph_traits<GraphSecond>::null_vertex()){
    			ita = prot1->tm_chains.find(prot1->get_chains()[vertex1]);
    			itb = prot2->tm_chains.find(prot2->get_chains()[get(correspondence_map_1_to_2, vertex1)]);
    			// Match both tm and globular chains
    			if((ita != prot1->tm_chains.end() && itb != prot2->tm_chains.end())||(ita == prot1->tm_chains.end() && itb == prot2->tm_chains.end())){
    				//cout << "Chains are both glob or both tm" << endl;
    				t++;
    				// Score for tm chains
    				if(ita != prot1->tm_chains.end() && itb != prot2->tm_chains.end()){
    					//tm-tm
						//cout << ita->second << " \t " << itb->second << endl;
						score += abs(ita->second-itb->second);
    				}else{
    					//cout << "glob-glob" << endl;
    					score += 5;
    				}

    			}else{
    				//cout << "mismatch" << endl;
    				score += 100;
    			}
    		}
    	}    	
    	//cout << "Score = " << score << endl;
    	//if(t > *q){
    		//*q = t;
    	if(score < *init_score){
    		//cout << "New score:\t" << score << " (was " << *init_score << ")\t" << "vertices: " << t << " (was " << *q << ")"<< endl;
    		*q = t;
    		*init_score = score;
	    	int i = 0;
	    	BGL_FORALL_VERTICES_T(vertex1, m_graph1, GraphFirst) {
	    		if (get(correspondence_map_1_to_2, vertex1) != graph_traits<GraphSecond>::null_vertex()){
	    			ita = prot1->tm_chains.find(prot1->get_chains()[vertex1]);
	    			itb = prot2->tm_chains.find(prot2->get_chains()[get(correspondence_map_1_to_2, vertex1)]);
					if((ita != prot1->tm_chains.end() && itb != prot2->tm_chains.end())||(ita == prot1->tm_chains.end() && itb == prot2->tm_chains.end())){
						correspond1[i] = vertex1;	
						correspond2[i] = get(correspondence_map_1_to_2, vertex1);
						i++;
	    			}
	    		}
	    	}
    	}		
    	//}
	    return (true);
  	}

private:
	int* q;
	int* n;
	int* init_score;
	int* correspond1;
	int* correspond2;	
	map<string,int>::const_iterator ita, itb;
	PDB* prot1;
	PDB* prot2;	
    const GraphFirst& m_graph1;
    const GraphSecond& m_graph2;
};

int mcs(PDB* prot1, PDB* prot2, unsigned int max_c){

	if(num_vertices(prot1->g) == 1){
		cout << endl << prot1->pdbname << " only contains one chain - skipping MCS search." << endl << endl;
		return(0);
	}
	if(num_vertices(prot2->g) == 1){
		cout << endl << prot2->pdbname << " only contains one chain - skipping MCS search." << endl << endl;
		return(0);
	}
	if(num_vertices(prot1->g)+num_vertices(prot2->g) > max_c){
		cout << endl << "More than " << max_c << " chains in total (" << num_vertices(prot1->g)+num_vertices(prot2->g) << ") - skipping MCS search." << endl << endl;
		return(0);
	}

	time_t start,end;
	int* p = new int(0);
	int* m = new int(0);
	int* s = new int(10000);
	int* c1 = new int[max(num_vertices(prot1->g),num_vertices(prot2->g))];
	int* c2 = new int[max(num_vertices(prot1->g),num_vertices(prot2->g))];
	map<string,int>::const_iterator it;
	
	cout << endl << "Searching for maximum common subgraph..." << endl;
	print_callback<Graph, Graph> my_callback(prot1, prot2, p, m, c1, c2, s);
	time (&start);
	
	//
	//mcgregor_common_subgraphs_maximum_unique(prot1->g, prot2->g, true, my_callback,vertices_equivalent(make_property_map_equivalent(prot1->vertexHomologyGroupMap,prot2->vertexHomologyGroupMap)));
	mcgregor_common_subgraphs_maximum_unique(prot1->g, prot2->g, true, my_callback);
	
	time (&end);
	double dif = difftime (end,start);
	cout << *m << " common subgraphs found in " << dif << " second";
	if(dif != 1) cout << "s";
	cout << ", best match contains " << *p << " vertices." << endl;
	//cout << "Best score: " << *s << endl;

	if(num_vertices(prot1->g) > num_vertices(prot2->g) && num_vertices(prot2->g) == (unsigned int)*p){
		cout << prot2->pdbname << " is contained within " << prot1->pdbname << endl;
	}
	if(num_vertices(prot2->g) > num_vertices(prot1->g) && num_vertices(prot1->g) == (unsigned int)*p){
		cout << prot1->pdbname << " is contained within " << prot2->pdbname << endl;
	}
	if(num_vertices(prot2->g) == num_vertices(prot1->g) && num_vertices(prot1->g) == (unsigned int)*p){
		cout << prot1->pdbname << " and " << prot2->pdbname << " are identical." << endl;
	}	

	if(*p){

		cout << endl << "Chain correspondences:" << endl;
		for(int i = 0; i < *p; i++){
		   cout << prot1->pdbname << ":" << prot1->get_chains()[c1[i]] << " ---> " << prot2->pdbname << ":" << prot2->get_chains()[c2[i]];
		   it = prot1->tm_chains.find(prot1->get_chains()[c1[i]]);
		   if(it != prot1->tm_chains.end()){
		   		cout << " (transmembrane/";
		   		// Add helix count match score here
		   }else{
		   		cout << " (globular/";
		   }
		   it = prot2->tm_chains.find(prot2->get_chains()[c2[i]]);
		   if(it != prot2->tm_chains.end()){
		   		cout << "transmembrane)" << endl;
				// Add helix count match score here
		   }else{
		   		cout << "globular)" << endl;
		   }
		   // Colour corresponding nodes
		   agsafeset(prot1->gv_nodes[c1[i]], (char*)"color", (char*)"orange", (char*)"");
		   agsafeset(prot2->gv_nodes[c2[i]], (char*)"color", (char*)"orange", (char*)"");	

		   prot1->correspond.push_back(c1[i]);
		   prot2->correspond.push_back(c2[i]);
		}

		Agnode_t *v, *h, *t;
		Agedge_t *e;

		for(v = agfstnode(prot1->gv); v; v = agnxtnode(prot1->gv,v)){
			for(e = agfstout(prot1->gv,v); e; e = agnxtout(prot1->gv,e)){
				h = aghead(e);
				t = agtail(e);
				char* tmp1 = agget(h,(char*)"name");
				char* tmp2 = agget(t,(char*)"name");
				bool cn1 = false, cn2 = false;
				for(int i = 0; i < *p; i++){
					if(prot1->get_chains()[c1[i]] == tmp1){
						cn1 = true;
					}	
					if(prot1->get_chains()[c1[i]] == tmp2){
						cn2 = true;
					}	
				}
				// Colour corresponding edges
				if(cn1 && cn2){
					agsafeset(e, (char*)"color", (char*)"orange", (char*)"");
				}	
			}		
		}

		for(v = agfstnode(prot2->gv); v; v = agnxtnode(prot2->gv,v)){
			for(e = agfstout(prot2->gv,v); e; e = agnxtout(prot2->gv,e)){
				h = aghead(e);
				t = agtail(e);
				char* tmp1 = agget(h,(char*)"name");
				char* tmp2 = agget(t,(char*)"name");
				bool cn1 = false, cn2 = false;
				for(int i = 0; i < *p; i++){
					if(prot2->get_chains()[c2[i]] == tmp1){
						cn1 = true;
					}	
					if(prot2->get_chains()[c2[i]] == tmp2){
						cn2 = true;
					}	
				}
				// Colour corresponding edges
				if(cn1 && cn2){
					agsafeset(e, (char*)"color", (char*)"orange", (char*)"");
				}	
			}		
		}
	}

	cout << endl;
	delete m;
	delete p;
	delete s;
	delete [] c1;
	delete [] c2;
	return(1);

}

void print_header(char* gvc_version){
	cout << endl << "Build and compare graphs generated from protein complexes" << endl << endl;
	cout << "(c) Timothy Nugent 2013, version " << VERSION << endl;
	cout << "Compiled using Boost " << BOOST_LIB_VERSION << " and Graphviz " << gvc_version << endl << endl;	
}

void usage(const char* progname, char* gvc_version){
	print_header(gvc_version);
 	cout << "Usage : " << progname << " [options] <PDB file 1> [PDB file 2]" << endl << endl;
 	cout << "Constructs a graph where chains are vertices and interchain contacts form" << endl;
 	cout << "edges. If two PDB files are provided, a maximum common subgraph algorithm" << endl;
 	cout << "will be run with the common subgraph annotated on the output graphs." << endl << endl;
	cout << "Options:" << endl;
	cout << "-a <string> Transmembrane topology of first PDB file, e.g. A:10,B:6,E:3" << endl;
	cout << "            A = transmembrane chain ID:10 = number of transmembrane helices, etc." << endl;
	cout << "-b <string> Transmembrane topology of second PDB file" << endl;
	cout << "-t <int>    Process transmembrane chains only. default 0" << endl;
	cout << "-l <int>    Layout algorithm. 1 = sfdp (default), 2 = fdp, 3 = circo, 4 = neato, 5 = kamada-kawai" << endl;
	//cout << "            Only circo and neato algorithms are deterministic" << endl;
	cout << "-d <int>    Draw graphs. 1 = both (default), 2 = PDB 1, 3 = PDB 2, 4 = none" << endl;
	cout << "-w <int>    Draw maximum common subgraphs. 1 = both (default), 2 = PDB 1, 3 = PDB 2, 4 = none" << endl;
	cout << "-m <int>    Maximum (combined) number of vertices/chains to run MCS algorithm - complexity is O(V1*V2). default 35" << endl;
	cout << "-i <float>  Sequence identity threshold for interchain homology match (uses Smith–Waterman). default " << SEQ_ID_THRESHOLD << endl;
	cout << "-s <int>    Perform symmetry search. 1 = slow (homology and subgraph match), 2 = fast (homology). default 1" << endl;
	cout << "-e <int>    Maximum number of vertices/chains to run symmetry search. default 27" << endl;
	cout << "-o <string> Output directory" << endl;
	cout << "-h          Display usage." << endl << endl;
}	

void split(string& arg, char_separator<char> sep, vector<string>& tmp){
    tokenizer< char_separator<char> > tokens(arg, sep);
    BOOST_FOREACH (const string& t, tokens){
    	tmp.push_back(t);
    }
}

void tokenize_chains(const char* arg, PDB* prot){
	
	char_separator<char> com(","), col(":");
	string c = arg;
	vector<string> tm_chains;
	split(c,com,tm_chains);
	BOOST_FOREACH (const string& t, tm_chains){
		c = t;
		vector<string> topology;
		split(c,col,topology);
		if(topology.size() == 2){
			//cout << topology[0] << " ---> " << topology[1] << endl;
			prot->set_tm_chains(topology[0],lexical_cast<int>(topology[1]));
		}
	}
	//cout << endl;	
}

void cleanup(PDB* prot1, PDB* prot2, GVC_t* gvc){

	delete prot1;
	delete prot2;
	gvFreeContext(gvc);
}

int main(int argc, const char* argv[]){

	PDB* prot1 = new PDB();
	PDB* prot2 = new PDB();
	GVC_t* gvc = gvContext();	

	int i = 1, w = 1, d = 1, m = 35, sym = 0, sym_max = 27, tm_only = 0;
	double s = SEQ_ID_THRESHOLD;
	if (argc < 2){
		usage(argv[0],gvcVersion(gvc));
		cleanup(prot1,prot2,gvc);
		return(1);
	}

	print_header(gvcVersion(gvc));

	while(i < argc){
  		if(argv[i][0] == '-'){
    		i++;
    		switch(argv[i-1][1]){
    			case 'a' : {tokenize_chains(argv[i],prot1); break;}
    			case 'b' : {tokenize_chains(argv[i],prot2); break;}
    			case 'o' : {prot1->output = argv[i];prot2->output = argv[i]; break;}
    			case 'l' : {prot1->layout = atoi(argv[i]);prot2->layout = atoi(argv[i]); break;}
    			case 'm' : {m = atoi(argv[i]); break;}
    			case 'w' : {w = atoi(argv[i]); break;}
    			case 'd' : {d = atoi(argv[i]); break;}
    			case 'i' : {s = atof(argv[i]); break;}
    			case 's' : {sym = atoi(argv[i]); break;}
    			case 'e' : {sym_max = atoi(argv[i]); break;}
    			case 't' : {tm_only = atoi(argv[i]); break;}
    			case 'h' : {usage(argv[0],gvcVersion(gvc));cleanup(prot1,prot2,gvc);return(1);}
    			default  : {usage(argv[0],gvcVersion(gvc));cleanup(prot1,prot2,gvc);return(1);}
			}   
   		}
   		i++;
   	}

   	cout << "Called with:" << endl;
   	for(int i = 0; i < argc; i++){
   		cout << argv[i] << " ";
   	}
   	cout << endl << endl;

   	boost::regex e("\\.pdb$",boost::regex::perl|boost::regex::icase);
   	if(regex_search(argv[i-1],e) && regex_search(argv[i-2],e)){

		prot1->parse_pdb(argv[i-2],s,tm_only);
		prot1->build_graph();  
		cout << endl; 		
		prot2->parse_pdb(argv[i-1],s,tm_only);
		prot2->build_graph();  
		int sg = mcs(prot1,prot2,m);	
		int l1 = prot1->layout_graph(gvc, true);	
		int l2 = prot2->layout_graph(gvc, false);
		if(sym) prot1->find_symmetry(sym,sym_max,1);
		if(sym) prot2->find_symmetry(sym,sym_max,1);
		if(l1 && (d == 1 || d == 2)) prot1->write_graph(gvc);	
		if(l2 && (d == 1 || d == 3)) prot2->write_graph(gvc);
		if(sg && (w == 1 || w == 2)) prot1->write_subgraph(gvc, prot2->pdbname);	
		if(sg && (w == 1 || w == 3)) prot2->write_subgraph(gvc, prot1->pdbname);
		if(l1) prot1->free_layout(gvc);
		if(l2) prot2->free_layout(gvc);	

   	}else if(regex_search(argv[i-1],e)){

		prot1->parse_pdb(argv[i-1],s,tm_only);
		prot1->build_graph();  
		cout << endl; 		
		int l1 = prot1->layout_graph(gvc, true);
		if(sym) prot1->find_symmetry(sym,sym_max,0);
		if(l1 && (d == 1 || d == 2)) prot1->write_graph(gvc);
		if(l1) prot1->free_layout(gvc);

   	}else{

   		cout << "Couldn't find any PDB files!" << endl;
   	
   	}
   	cout << endl;

   	cleanup(prot1,prot2,gvc);

	return(0);
}
