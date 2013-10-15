#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

namespace alignment {

	double smithwaterman(string& seq_a, string& seq_b, double match, double mismatch, double indel, bool verbose){

		int m = seq_a.length();
	  	int n = seq_b.length();

		// Initialise matrix H
		double** H = new double*[m+1];
		for(int i = 0; i < m+1; i++){
			H[i] = new double[n+1];
			for(int j = 0; j <= n; j++){
				H[i][j] = 0.0;	
			}	
		}

		// Initialise matrices to store path
		int** Trace_a = new int*[m+1];
		int** Trace_b = new int*[m+1];
		for(int i = 0; i < m+1; i++){
		   Trace_a[i] = new int[n+1];
		   Trace_b[i] = new int[n+1];
		}

		int idx = 0;
		double temp[4];
		for(int i=1; i <= m; i++){
			for(int j = 1; j <= n; j++){

				// http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
				if(seq_a[i-1] == seq_b[j-1]){
					temp[0] = H[i-1][j-1] + match;
				}else{
					temp[0] = H[i-1][j-1] + mismatch;
				}
				temp[1] = H[i-1][j] + indel;                  
				temp[2] = H[i][j-1] + indel;                 
				temp[3] = 0.0;

				idx = 0;
			  	for(int l = 0; l < 4; l++){
					if(temp[l] > H[i][j]){
						H[i][j] = temp[l];
						idx = l; 
					}
			  	}

				switch(idx){
		      		case 0:  
		      			// Match or mismatch
		   				Trace_a[i][j] = i-1;
						Trace_b[i][j] = j-1;
						break;
		      		case 1:     
		      		    // Deletion in seq_a
		     			Trace_a[i][j] = i-1;
						Trace_b[i][j] = j;
						break;
		      		case 2:         
		      			// Deletion in seq_b
		      			Trace_a[i][j] = i;
						Trace_b[i][j] = j-1;
						break;
		      		case 3:      
		      			// Subsequence
		      			Trace_a[i][j] = i;
						Trace_b[i][j] = j;	
					break;
		      	}
		    }
		}

		// Search H for the maximum score
		double H_max = 0.;
	  	int curr_i = 0, curr_j = 0;
	  	for(int i = 1; i <= m; i++){
	    	for(int j = 1; j <= n; j++){
	      		if(H[i][j] > H_max){
					H_max = H[i][j];
					curr_i = i;
					curr_j = j;
	      		}
	    	}
	  	}
	    
	    int pos = 0, ident = 0;
		int next_i = Trace_a[curr_i][curr_j];
		int next_j = Trace_b[curr_i][curr_j];
		char *aligned_a = new char[m+n+2];
		char *aligned_b = new char[m+n+2];

		// Traceback from H_max
		while(((curr_i != next_i) || (curr_j != next_j)) && (next_j != 0) && (next_i != 0)){

		    if(next_i==curr_i){
		    	// gap in A
		    	aligned_a[pos] = '-';                  
		    }else{
		    	// match/mismatch in A
		    	aligned_a[pos] = seq_a[curr_i-1];   
		    }                   

		    if(next_j==curr_j){
		    	// gap in B
		    	aligned_b[pos] = '-';                  
		    }else{	
		    	// match/mismatch in B
		    	aligned_b[pos] = seq_b[curr_j-1]; 
		    	if(aligned_b[pos] == aligned_a[pos]){
		    		ident++;	
		    	} 
		    }

		    curr_i = next_i;
		    curr_j = next_j;
		    next_i = Trace_a[curr_i][curr_j];
		    next_j = Trace_b[curr_i][curr_j];
		    pos++;
	    }
	  	
	  	if(verbose){
	  		cout << endl;
		  	for(int i = pos-1; i >= 0; i--) cout << aligned_a[i]; 
		  	cout << endl;
		  	for(int j = pos-1; j >= 0; j--) cout << aligned_b[j];
		  	cout << endl;
		}

		for(int i = 0; i < m+1; ++i){
		   delete [] Trace_a[i];
		   delete [] Trace_b[i];
		   delete [] H[i];
		}
		delete [] Trace_a;
		delete [] Trace_b;
		delete [] H;
		delete [] aligned_a;
		delete [] aligned_b;

		// Return sequence identity
		if(ident){
	  		return((double)ident/(min(seq_a.length(),seq_b.length())));
	  	}else{
	  		return(0.0);
	  	}	
	}
}
