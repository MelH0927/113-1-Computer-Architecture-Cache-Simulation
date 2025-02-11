#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <iterator>

using namespace std;

class cache_tag
{
	private:
		vector<bool> tag;
		bool valid;
		int tag_size;
	
	public:
		cache_tag(int size)
		{
			tag_size = size;
			for (int i=0; i<size; i++)
				tag.push_back(0);
			valid = 0;
		}
		
		void setup_cache_tag(const vector<bool>& tag_ref)
		{
			valid = 1;
			for (int i=0; i<tag_size; i++)
				tag[i] = tag_ref[i];
		}
		
		bool cache_tag_hit(const vector<bool>& tag_ref)
		{
			if (valid == 0)
				return 0;
			for (int i=0; i<tag_size; i++)
				if (tag[i] != tag_ref[i])
					return 0;
			return 1;
		}
};

class cache_set
{
	private:
		vector<cache_tag> set_array;
		set<int> least_used;
		int associativity;

		void refill_least_used(int item)
		{
			if (associativity > 1 && least_used.size() == 1 && *least_used.begin() == item)
				for (int j=0; j<associativity; j++)
					least_used.insert(j);
		}
	
	public:
		cache_set(int asso, int tag_size)
		{
			cache_tag tag_instance(tag_size);
			for (int i=0; i<asso; i++) {
				set_array.push_back(tag_instance);
				least_used.insert(i);
			}
			associativity = asso;
		}

		bool cache_set_hit(const vector<bool>& tag_ref)
		{
			// hit
			for (int i=0; i<associativity; i++) {
				if (set_array[i].cache_tag_hit(tag_ref)) {
					refill_least_used(i);
					if (associativity > 1 && least_used.count(i))
						least_used.erase(i);
					return 1;
				}
			}
			
			// miss
			int pick = *least_used.begin();
			refill_least_used(pick);
			set_array[pick].setup_cache_tag(tag_ref);
			if (associativity > 1 && least_used.count(pick))
				least_used.erase(pick);
			return 0;
		}
};

class cache
{
	private:
		vector<cache_set> cache_sets;
		set<int> indexing_bits;
		int address_size, set_count, index_size, tag_size;
	
	public:
		cache(int Address_bits, int Block_size, int Cache_sets, int Associativity, 
			const set<int>& idx_bits_ref, ofstream& outfile)
		{
			// data setup
			int offset_bit_count = log2(Block_size);
			index_size = log2(Cache_sets);
			set<int>::iterator it=idx_bits_ref.begin();
			for (; it != idx_bits_ref.end(); it++)
				indexing_bits.insert(*it);
			address_size = Address_bits;
			set_count = Cache_sets;
			tag_size = address_size - index_size - offset_bit_count;
			
			// cache storage setup
			cache_set set_instance(Associativity, tag_size);
			for (int i=0; i<set_count; i++)
				cache_sets.push_back(set_instance);
			
			// print info
			outfile << "Address bits: " << Address_bits << endl;
			outfile << "Block size: " << Block_size << endl;
			outfile << "Cache sets: " << Cache_sets << endl;
			outfile << "Associativity: " << Associativity << endl;
			outfile << "\n\n";
			outfile << "Offset bit count: " << offset_bit_count << endl;
			outfile << "Indexing bit count: " << index_size << endl;
			outfile << "Indexing bits:";
			it=indexing_bits.begin();
			for (; it != indexing_bits.end(); it++) {
				outfile << " " << Address_bits - *it - 1;
				//cout << *it << endl;
			}
			outfile << "\n\n";
		}

		bool cache_hit(const vector<bool>& addr)
		{
			int index = 0;
			vector<bool> tag;
			for (int i=0; i<(index_size + tag_size); i++){
				if (indexing_bits.count(i)) {
					index *= 2;
					index += addr[i];
				}
				else
					tag.push_back(addr[i]);
			}
			return cache_sets[index].cache_set_hit(tag);
		}
};

void find_index_bits(set<int>& indexing_bits, const vector<vector<bool> >& ref,
			int address_bits, int ref_counts, int idx_count, int offset)
{
	// lsb
	/*for (int i=0; i<idx_count; i++)
		indexing_bits.insert(i+offset);*/
	
	// zero cost
	
	// setup corr matrix
	vector <vector<double> > corr_matrix;
	vector <double> corr_row;
	int idxtag_size = address_bits - offset;
	for (int i=0; i<idxtag_size; i++)
		corr_row.push_back(0);
	for (int i=0; i<idxtag_size; i++)
		corr_matrix.push_back(corr_row);
	
	// count E_i, j
	for (int i=0; i<ref_counts; i++)
		for (int j=0; j<idxtag_size; j++)
			for (int k=j+1; k<idxtag_size; k++)
				if (ref[i][j] == ref[i][k]) {
					corr_matrix[j][k] += 1;
					corr_matrix[k][j] += 1;
				}
	
	// cal C_i, j
	for (int i=0; i<idxtag_size; i++) {
		for (int j=0; j<idxtag_size; j++) {
			int E = corr_matrix[i][j];
			int D = ref_counts - E;
			if (E > D)
				corr_matrix[i][j] = (double) D/E;
			else
				corr_matrix[i][j] = (double) E/D;
		}
	}
	
	
	// setup quality meas
	vector <double> qual_meas;
	for (int i=0; i<idxtag_size; i++) {
		int Z = 1;
		int O = 1;
		for (int j=0; j<ref_counts; j++)
			if (ref[j][i] == 0)
				Z += 1;
			else
				O += 1;
		if (Z<O)
			qual_meas.push_back((double) Z/O);
		else
			qual_meas.push_back((double) O/Z);
	}

	// get indexing bit
	for (int i=0; i<idx_count; i++) {
		double Q = 0;
		int idx = idxtag_size-1;
		while (indexing_bits.count(idx))
			idx--;
		for (int j=0; j<idxtag_size; j++) {
			//cout << j << ": ";
			if (qual_meas[j] > Q) {
				//cout << "y";
				Q = qual_meas[j];
				idx = j;
			}
			//cout << endl;
		}
		indexing_bits.insert(idx);
		qual_meas[idx] = 0;
		for (int j=0; j<idxtag_size; j++) {
			qual_meas[j] *= corr_matrix[j][idx];
		}
		//cout << endl;
		//cout << idx << endl;
	}
	/*set<int>::iterator it = indexing_bits.begin();
	for (; it != indexing_bits.end(); it++)
		cout << *it << endl;*/
}

int main(int argc, char *argv[])
{
	ifstream infile;
	string info;
	int Address_bits, Block_size, Cache_sets, Associativity;

	ofstream outfile;
	string line1, line2, data;
	vector<vector<bool> > ref;
	vector<bool> ref_instance;
	
	set<int> indexing_bits;
	
	int miss = 0;
	string result;

	// read config
	infile.open(argv[1], ios::in);
	infile >> info >> Address_bits;
	infile >> info >> Block_size;
	infile >> info >> Cache_sets;
	infile >> info >> Associativity;
	infile.close();

	// read testbench
	for (int i=0; i<Address_bits; i++)
		ref_instance.push_back(0);
	infile.open(argv[2], ios::in);
	getline(infile, line1);
	getline(infile, data);
	while (data != ".end") {
		for (int i=0; i<Address_bits; i++)
			ref_instance[i] = (data[i] == '1');
		ref.push_back(ref_instance);
		getline(infile, data);
	}
	line2 = data;
	infile.close();
	
	// find index bits
	find_index_bits(indexing_bits, ref, Address_bits, ref.size(), log2(Cache_sets), log2(Block_size));
	
	// setup cache
	outfile.open(argv[3], ios::out);
	cache mycache(Address_bits, Block_size, Cache_sets, Associativity, indexing_bits, outfile);

	// find hit or miss
	outfile << line1 << endl;
	for (int i=0; i<ref.size(); i++) {
		for (int j=0; j<Address_bits; j++)
			outfile << ref[i][j];
		if (mycache.cache_hit(ref[i]))
			result = "hit";
		else {
			result = "miss";
			miss += 1;
		}
		outfile << " " << result << endl;
	}
	outfile << line2 << "\n\nTotal cache miss count: " << miss << endl;
	outfile.close();

	return 0;
}
