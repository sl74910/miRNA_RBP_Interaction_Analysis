#include <bits/stdc++.h>

using namespace std;
#define int long long

// Parse each line in the file and return a vector containing each field
vector<string> parse_line(const string& line) {
    stringstream ss(line);
    string field;
    vector<string> fields;
    while (ss >> field) {
        fields.push_back(field);
    }
    return fields;
}

// Merge RBP and miRNA mRNA interaction information, remove rows where RBP and mRNA genes are the same, and output two intermediate files
void merge_mrna_interactions(const string& rbp_file, const string& mirna_file, const string& output_file1, const string& output_file2) {
    unordered_set<string> rbp_genes;  // Used to store RBP gene information
    ifstream in_rbp(rbp_file);
    ifstream in_mirna(mirna_file);
    ofstream out_rbp(output_file1);
    ofstream out_mirna(output_file2);

    string line;

    // Check if the files opened successfully
    if (!in_rbp || !in_mirna) {
        cout << "Error: Unable to open input files." << endl;
        return;
    }

    // Process the RBP file, filter out rows where RBP gene and target gene are the same, and write valid rows to the output file
    while (getline(in_rbp, line)) {
        vector<string> fields = parse_line(line);
        if (fields.size() < 7) continue;  // Skip incomplete rows

        string s4 = fields[3];  // RBP gene
        string s7 = fields[6];  // Target gene

        // Skip rows where RBP gene is the same as the target gene
        if (s4 == s7) continue;

        rbp_genes.insert(s4);  // Store RBP gene in the set
        out_rbp << line << endl;  // Write valid RBP rows to the output file
    }

    // Process the miRNA file and filter out miRNA rows related to RBP
    int i = 0;
    while (getline(in_mirna, line)) {
        vector<string> fields = parse_line(line);
        if (fields.size() < 7) continue;  // Skip incomplete rows

        string s7 = fields[6];  // miRNA target gene

        // If the gene does not appear in the RBP file, write it to the output file
        if (rbp_genes.find(s7) == rbp_genes.end()) {
            out_mirna << line << endl;
        }

        i++;
        if (i % 100 == 0) cout << "Processed " << i << " miRNA entries." << endl;
    }

    cout << "Merging and filtering completed." << endl;
}

// Merge RBP and miRNA mRNA interaction information
void merge_rbp_mirna_interactions(const string& rbp_file, const string& mirna_file, const string& output_file) {
    map<string, vector<string>> rbp_map;  // Used to store RBP-mRNA interaction information
    ifstream in_rbp(rbp_file);
    ifstream in_mirna(mirna_file);
    ofstream out(output_file);

    string line;

    // Check if the files opened successfully
    if (!in_rbp || !in_mirna) {
        cout << "Error: Unable to open input files." << endl;
        return;
    }

    // Read the RBP file and store it in the map
    while (getline(in_rbp, line)) {
        stringstream ss(line);
        string extra_info, start, end, rbp, xa1, xa2, mRNA;
        ss >> extra_info >> start >> end >> rbp >> xa1 >> xa2 >> mRNA;
        string combined_info = rbp + '\t' + mRNA + '\t' + extra_info + '\t' + start + '\t' + end;
        rbp_map[mRNA].push_back(combined_info);  // Store by mRNA classification
    }

    cout << "Step 1: Merging miRNA-mRNA and RBP-mRNA interactions..." << endl;

    int line_count = 0;
    // Read the miRNA file and check for corresponding RBP interactions
    while (getline(in_mirna, line)) {
        stringstream ss(line);
        string extra_info, start, end, miRNA, xa1, xa2, mRNA;
        ss >> extra_info >> start >> end >> miRNA >> xa1 >> xa2 >> mRNA;
        string combined_info = miRNA + '\t' + mRNA + '\t' + extra_info + '\t' + start + '\t' + end;
        // If this miRNA's target mRNA has RBP interactions
        if (rbp_map.count(mRNA)) {
            for (const auto& rbp_info : rbp_map[mRNA]) {
                out << combined_info << '\t' << rbp_info << endl;
            }
        }
        
        line_count++;
        if (line_count % 100 == 0) {
            cout << "Processed " << line_count << " miRNA-mRNA pairs." << endl;
        }
    }

    cout << "Merging completed." << endl;
}

// Analyze the overlap and adjacency between RBP and miRNA interactions
void analyze_interactions(const string& input_file, const string& output_file) {
    ifstream in(input_file);
    ofstream out(output_file);
    string line;

    // Maps to count information about total, overlapping, and adjacent interactions
    map<string, int> total_count, overlap_count, adjacent_count;

    if (!in) {
        cout << "Error: Unable to open the input file." << endl;
        return;
    }

    int line_count = 0;
    while (getline(in, line)) {
        stringstream ss(line);
        string miRNA, mRNA, extra_info1, rbp, mRNA2, extra_info2;
        int start_rbp, end_rbp, start_mirna, end_mirna;
        ss >> miRNA >> mRNA >> extra_info1 >> start_mirna >> end_mirna >> rbp >> mRNA >> extra_info2 >> start_rbp >> end_rbp;
        string interaction_key = rbp + '\t' + mRNA + '\t' + miRNA;  // Unique RBP-miRNA combination
        total_count[interaction_key]++;  // Count total interactions

        // Determine whether the interactions overlap or are adjacent
        if (extra_info1 == extra_info2 && 
            abs(end_rbp - start_rbp) + 1 + abs(end_mirna - start_mirna) + 1 > 
            max({abs(end_rbp - start_mirna), abs(end_rbp - end_mirna), abs(start_rbp - start_mirna), abs(start_rbp - end_mirna)}) + 1 ) {
            overlap_count[interaction_key]++;
        } else {
            adjacent_count[interaction_key]++;
        }

        line_count++;
        if (line_count % 100000 == 0) {
            cout << "Processed " << line_count << " interactions." << endl;
        }
    }

    // Output analysis results
    for (const auto& entry : total_count) {
        out << entry.first << '\t' << overlap_count[entry.first] << '\t' << adjacent_count[entry.first] << endl;
    }

    cout << "Analysis completed." << endl;
}

signed main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    string rbp_file = "./inputs/RBPmRNAv27.bed";
    string mirna_file = "./inputs/miRNAmRNAv27.bed";
    string intermediate_file1 = "./outputs/RBPmRNAv27_filtered.bed";
    string intermediate_file2 = "./outputs/miRNAmRNAv27_filtered.bed";
    string merged_file = "./outputs/RBP_m_mi_intermediate.bed";
    string output_file = "./outputs/RBP_m_mi.bed";

    // Step 1: Filter RBP and miRNA data based on mRNA genes
    cout << "Step 1: Merging and filtering RBP and miRNA based on mRNA gene..." << endl;
    merge_mrna_interactions(rbp_file, mirna_file, intermediate_file1, intermediate_file2);
    
    // Step 2: Merge RBP and miRNA interaction information
    cout << "Step 2: Merging filtered RBP and miRNA interactions..." << endl;
    merge_rbp_mirna_interactions(intermediate_file1, intermediate_file2, merged_file);

    // Step 3: Analyze the adjacency and overlap between RBP and miRNA interactions
    cout << "Step 3: Analyzing adjacent and overlapping interactions..." << endl;
    analyze_interactions(merged_file, output_file);

    return 0;
}
