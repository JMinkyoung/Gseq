#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <vector>
#include <ctime>

using namespace std;

#define REF_LEN 50000
#define REF_FILE "referenceDNA_50000.txt"
#define K 70
#define N 100000
#define SNP_NUM 2
#define SHORT_READ_FILE "short_read_dna_change.txt"
#define SHORT_READ_ORIGIN_FILE "short_reads_seq_70_100million2.txt"

const char ALPHABET[4] = { 'a', 't', 'g', 'c' };

// 모든 target dna sequence를 빠짐 없이 나타낼 수 있는 n개의 short_read를 생성하고 파일로 저장한다.
// n개의 short read로 모든 dna sequence를 표현하려면 최소 dna_sequence 길이에 n을 나눈 길이만큼 순차적으로 샘플링해야한다.
// 따라서 k는 dna length / n보다 커야하고, dna length / n ~ k 사이의 간격만큼 랜덤으로 인덱스를 변경하여 short read를 샘플링하면 모든 sequence를 표현할 수 있다.
// 단 예외가 있는데 시작 sequence와 끝 sequence는 반드시 포함해야한다.
// 첫번째 iteration에는 모든 sequence를 포함하는 short read를 순차적으로 생성하고, 남은 횟수가 있으면 그 만큼 랜덤으로 샘플링한다.
void generate_short_read(string filename, string target_dna, int k, int n)
{
	// n개의 샘플링으로 모든 sequence를 표현하기 위한 short read의 최소 길이
	int minimum_hop = REF_LEN / n;
	if (minimum_hop > k)
	{
		cout << n << "개의 short read로는 모든 sequence를 포함할 수 없습니다." << endl;
		exit(0);
	}
	// REF_LEN / n이 나누어 떨어지지 않으면 minimum_hop에 1을 더해줌
	if (REF_LEN % n != 0)
		minimum_hop++;

	// minimum_hop~k 사이의 간격으로 샘플링하면 된다.
	int hop_range = k - minimum_hop;

	ofstream write_file;
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0toRange(0, hop_range);
	uniform_int_distribution<mt19937::result_type> from0toBefore_k(0, REF_LEN - k - 1);

	write_file.open(filename);


	int short_cnt = 1;
	int str_idx = 0;
	string short_read = target_dna.substr(REF_LEN - k, k);
	write_file << short_read << endl;
	// 모든 sequence가 포함되게 short_read를 생성. 겹치는 부분이 랜덤적으로 존재 (범위 rand_hop)
	while (short_cnt < n)
	{
		short_read = target_dna.substr(str_idx, k);
		int rand_hop = from0toRange(mersenne);
		str_idx += rand_hop + minimum_hop;
		short_cnt++;
		write_file << short_read;
		// 더 이상 샘플링 할 수 없으면 break함.
		if (str_idx + k > REF_LEN)
			break;
		write_file << endl;
	}
	// 남은 short_cnt에 대해 랜덤으로 위치를 골라 k개의 sequence 생성
	while (short_cnt < n)
	{
		write_file << endl;
		str_idx = from0toBefore_k(mersenne);
		short_read = target_dna.substr(str_idx, k);
		short_cnt++;
		write_file << short_read;
	}

	write_file.close();

}

// reference DNA에서 선형적으로 k길이 만큼 2개 정도 SNP가 있다고 하였을때 다른 DNA sequence를 만들어준다.
// SNP의 위치는 0~k사이의 길이 중 2개가 있으나 위치는 랜덤이며, reference DNA과 동일한 염기를 가질 수 있기 때문에 k개 sequence에 대해 0~4개 정도의 SNP가 존재하게 된다.
string change_snp_reference_sequence(string ref_dna, int k)
{
	// 이미 생성한 파일이 있으면 해당 파일을 읽어옴
	ifstream read_file;
	read_file.open(SHORT_READ_ORIGIN_FILE);
	string return_str = "";
	if (read_file.is_open())
	{
		while (!read_file.eof())
		{
			string str;
			getline(read_file, str);
			return_str += str;
		}
		read_file.close();
		return return_str;
	}

	int sub_k_pos = 0;
	int dna_size = ref_dna.size();
	string dif_dna = string(ref_dna);
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0tok(0, k - 1);
	uniform_int_distribution<mt19937::result_type> from0to3(0, 3);

	int change_pos1, change_pos2;

	while (sub_k_pos + k < dna_size)
	{
		change_pos1 = from0tok(mersenne);
		change_pos2 = from0tok(mersenne);
		// 동일한 위치가 선택되었을 때 다른 위치를 고를때까지 다시 위치 선택
		while (change_pos1 == change_pos2)
			change_pos2 = from0tok(mersenne);
		// SNP 위치가 랜덤으로 결정되었으므로 값을 바꿔준다. 
		dif_dna[change_pos1 + sub_k_pos] = ALPHABET[from0to3(mersenne)];
		dif_dna[change_pos2 + sub_k_pos] = ALPHABET[from0to3(mersenne)];
		sub_k_pos += k;
	}
	int remains = dna_size - sub_k_pos;
	int snp_ratio = k / SNP_NUM;

	// DNA의 마지막의 애매하게 남은 부분에 대해 SNP를 변경
	for (int i = 0; i < remains / snp_ratio + 1; i++)
	{
		change_pos1 = from0tok(mersenne) % remains;
		dif_dna[change_pos1 + sub_k_pos] = ALPHABET[from0to3(mersenne)];
	}
	// short read를 합한 sequence와 비교하기 위해 파일로 저장
	ofstream write_file;
	write_file.open(SHORT_READ_ORIGIN_FILE);
	write_file.write(dif_dna.c_str(), dif_dna.size());
	write_file.close();
	return dif_dna;
}

string load_dna(string filename)
{
	// 이미 생성한 파일이 있으면 해당 파일을 읽어옴
	ifstream read_file;
	read_file.open(filename);
	string return_str = "";
	if (read_file.is_open())
	{
		while (!read_file.eof())
		{
			string str;
			getline(read_file, str);
			return_str += str;
		}
		read_file.close();
		return return_str;
	}
}

// 랜덤 염기서열을 생성하는 함수
string gen_random(string s, const int len)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dist(0, 3);

	static const char sequence[] = "atgc";

	// 새로운 난수 생성 방식
	for (int i = 0; i < len; i++) {
		s += sequence[dist(gen)];
	}
	return s;
}

// Reference DNA Sequence 파일 생성
void make_referenceDNA(int len) {
	ofstream fout;
	string str;

	fout.open("referenceDNA.txt");
	str = gen_random(str, 50000);
	fout << str;
	fout.close();
}
int main()
{

	// reference DNA 파일을 만든다.
	make_referenceDNA(50000);

	// my DNA와 Short Read 파일을 만든다.
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0to3(0, 3);
	string ref_dna = load_dna(REF_FILE);

	string target_dna = change_snp_reference_sequence(ref_dna, K);
	generate_short_read(SHORT_READ_FILE, target_dna, K, N);
}
