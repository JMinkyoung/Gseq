#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <vector>
#include <ctime>

using namespace std;

#define REF_LEN 10000
#define RULE_1 47
#define RULE_2 61
#define REF_FILE "referenceDNA.txt"
#define K 70
#define N 100000
#define SNP_NUM 2
#define SHORT_READ_FILE "short_read_dna_change.txt"
#define SHORT_READ_ORIGIN_FILE "recovered_dna.txt"


const char ALPHABET[4] = { 'A', 'T', 'G', 'C' };

// DNA sequence를 만든다. 문자열(DNA)의 반복을 줄이기 위해 47 배수 위치의 DNA 염기의 값을 이전 값과 다른 값을 넣어주고,
// 61 배수의 위치에는 DNA 염기 값을 A->T->G->C로 순서대로 넣어준다.
// 47과 61은 각각 소수로 서로 다른 변환과정이 겹치는 일이 없도록 하였고, 47*61의 배수의 위치의 경우 첫번째 룰을 우선시 하였다.
string generate_reference_dna(string filename, int size, mt19937 mersenne, uniform_int_distribution<mt19937::result_type> generator)
{
	string return_str = "";

	// 이미 생성한 파일이 있으면 해당 파일을 읽어옴
	ifstream read_file;
	read_file.open(filename);
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

	// 생성한 파일이 없으면 처음부터 생성
	int rule_2_pos = 0;
	return_str += ALPHABET[generator(mersenne)];
	for (int i = 1; i < size; i++)
	{
		// Rule 1 -> RULE_1 번째 위치마다 이전 위치의 염기의 반복을 없앰
		if (i % RULE_1 == 0)
		{
			int gen_idx = generator(mersenne);
			char last_char = return_str.back();
			while (last_char == ALPHABET[gen_idx])
			{
				gen_idx = generator(mersenne);
			}
			return_str += ALPHABET[gen_idx];
		}
		// Rule 2 -> RULE_2 번째 위치마다 미리 정해둔 순서대로 A->T->G->C 염기 할당
		else if (i % RULE_2 == 0)
		{
			return_str += ALPHABET[rule_2_pos];
			rule_2_pos = (rule_2_pos + 1) % 4;
		}
		else
			return_str += ALPHABET[generator(mersenne)];
	}
	ofstream write_file;
	write_file.open(filename);
	write_file.write(return_str.c_str(), return_str.size());
	write_file.close();
	return return_str;
}


// 모든 target dna sequence를 빠짐 없이 나타낼 수 있는 n개의 short_read를 생성하고 파일로 저장한다.
// n개의 short read로 모든 dna sequence를 표현하려면 최소 dna_sequence 길이에 n을 나눈 길이만큼 순차적으로 샘플링해야한다.
// 따라서 k는 dna length / n보다 커야하고, dna length / n ~ k 사이의 간격만큼 랜덤으로 인덱스를 변경하여 short read를 샘플링하면 모든 sequence를 표현할 수 있다.
// 단 예외가 있는데 시작 sequence와 끝 sequence는 반드시 포함해야한다.
// 첫번째 iteration에는 모든 sequence를 포함하는 short read를 순차적으로 생성하고, 남은 횟수가 있으면 그 만큼 랜덤으로 샘플링한다.
vector<string> generate_short_read(string filename, string target_dna, int k, int n)
{
	// 이미 생성한 파일이 있으면 해당 파일을 읽어옴
	ifstream read_file;
	vector<string> return_vec = vector<string>();
	read_file.open(filename);
	if (read_file.is_open())
	{
		while (!read_file.eof())
		{
			string str;
			getline(read_file, str);
			return_vec.push_back(str);
		}
		read_file.close();
		return return_vec;
	}

	// n개의 샘플링으로 모든 sequence를 표현하기 위한 short read의 최소 길이
	int minimum_hop = target_dna.size() / n;
	if (minimum_hop > k)
	{
		cout << n << "개의 short read로는 모든 sequence를 포함할 수 없습니다." << endl;
		exit(0);
	}
	// REF_LEN / n이 나누어 떨어지지 않으면 minimum_hop에 1을 더해줌
	if (target_dna.size() % n != 0)
		minimum_hop++;
	
	// minimum_hop~k 사이의 간격으로 샘플링하면 된다.
	int hop_range = k - minimum_hop;

	ofstream write_file;
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0toRange(0, hop_range);
	uniform_int_distribution<mt19937::result_type> from0toBefore_k(0, target_dna.size() -k-1);

	write_file.open(filename);
	
	
	int short_cnt = 1;
	int str_idx = 0;
	string short_read = target_dna.substr(target_dna.size() - k, k);
	return_vec.push_back(short_read);
	write_file << short_read << endl;
	// 모든 sequence가 포함되게 short_read를 생성. 겹치는 부분이 랜덤적으로 존재 (범위 rand_hop)
	while (short_cnt < n)
	{
		short_read = target_dna.substr(str_idx, k);
		int rand_hop = from0toRange(mersenne);
		str_idx += rand_hop + minimum_hop;
		short_cnt++;
		write_file << short_read;
		return_vec.push_back(short_read);
		// 더 이상 샘플링 할 수 없으면 break함.
		if (str_idx + k > target_dna.size())
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
		return_vec.push_back(short_read);
	}

	write_file.close();
	return return_vec;
	
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
	uniform_int_distribution<mt19937::result_type> from0tok(0, k-1);
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


// trivial method로 short read를 하나하나 비교해가며 threshold보다 mismatch가 작을 때 해당 short read의 위치를 구한다.
// mismatch threshold는 SNP를 변경한 dna sequence에서 k개의 길이 안에 최대 4개의 mismatch가 존재하도록 만들었기 때문에 4로 설정하였다.
string trivial_method(vector<string> short_reads, string ref_dna, int threshold, int k, int n)
{
	int mismatches = 0;
	int pos = 0;
	// DNA 염기서열을 복구하기 위해 reference dna를 복사한다.
	string recovered_dna = string(ref_dna);
	bool found = false;
	for (int i = 0; i < short_reads.size(); i++)
	{
		// debug용. short read의 위치를 찾지 못하였을 때 출력하고 프로그램을 종료함. 
		if (!found && i != 0)
		{
			cout << "cannot found " << i << endl;
		}
		found = false;
		// 하나의 short read를 전체 ref_dna에 대하여 k길이만큼의 sequence를 하나하나 position을 옮겨다니며 비교
		for (int j = 0; j < ref_dna.size() - k + 1; j++)
		{
			mismatches = 0;
			for (int l = 0; l < k; l++)
			{
				// mismatch가 있으면 1증가
				if (ref_dna[j + l] != short_reads[i][l])
					mismatches++;
				// mismatch가 threshold보다 크면 다음 위치로 이동
				if (mismatches > threshold)
					break;
			}
			// threshold보다 작은 mismatch가 나왔을 때 해당 위치를 short read의 dna위치로 보고 복구한다.
			if (mismatches <= threshold)
			{
				for (int l = 0; l < k; l++)
				{
					recovered_dna[j + l] = short_reads[i][l];
				}
				found = true;
				//cout << "found " << i << " at " << j << endl;
				break;
			}
		}
		if (i % 1000 == 0)
			cout << "iter = " << i << endl;
	}
	return recovered_dna;
}

int main()
{
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0to3(0, 3);
	string ref_dna = generate_reference_dna(REF_FILE, REF_LEN, mersenne, from0to3);
	
	string target_dna = change_snp_reference_sequence(ref_dna, K);
	vector<string> short_reads = generate_short_read(SHORT_READ_FILE, target_dna, K, N);
	clock_t start, end;
	double result;

	cout << "ref_dna size " << ref_dna.size() << endl;
	cout << "short read size " << short_reads.size() << endl;
	start = clock();
	string recovered_dna = trivial_method(short_reads, ref_dna, 2, K, N);
	end = clock();
	result = double(end - start);
	cout << "trivial method running time : " << result / CLOCKS_PER_SEC << "\n" << endl;

	ofstream write_file;
	write_file.open("recovered_dna.txt");
	write_file.write(recovered_dna.c_str(), recovered_dna.size());
	write_file.close();

	cout << "recovered_dna size " << recovered_dna.size() << endl;
	int cnt = 0;
	for (int i = 0; i < recovered_dna.size(); i++)
	{
		if (recovered_dna[i] == ref_dna[i])
			cnt++;
	}
	cout << "accuracy : " << (double)(1.0 * cnt / recovered_dna.size()) << endl;
}
