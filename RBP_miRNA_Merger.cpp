#include <bits/stdc++.h>

using namespace std;

// Parse each line of the file and return a vector containing the fields
vector<string> parse_line(const string& line) {
    stringstream ss(line);
    string field;
    vector<string> fields;
    while (ss >> field) {
        fields.push_back(field);
    }
    return fields;
}

// Merge RBP and miRNA mRNA interaction data, excluding rows where RBP and mRNA genes are the same, and output two intermediate files
void merge_mrna_interactions(const string& rbp_file, const string& mirna_file, const string& output_file1, const string& output_file2) {
    unordered_set<string> rbp_genes;  // To store RBP gene information
    ifstream in_rbp(rbp_file);
    ifstream in_mirna(mirna_file);
    ofstream out_rbp(output_file1);
    ofstream out_mirna(output_file2);

    string line;

    // Check if the input files were successfully opened
    if (!in_rbp || !in_mirna) {
        cerr << "Error: Unable to open input files." << endl;
        return;
    }

    // Process the RBP file, filter out rows where RBP and target genes are the same, and write the valid rows to the output file
    while (getline(in_rbp, line)) {
        vector<string> fields = parse_line(line);
        if (fields.size() < 7) continue;  // Skip incomplete rows

        string s4 = fields[3];  // RBP gene
        string s7 = fields[6];  // Target gene

        // Skip rows where RBP gene is the same as the target gene
        if (s4 == s7) continue;

        rbp_genes.insert(s4);  // Store the RBP gene in the set
        out_rbp << line << endl;  // Write the valid RBP row to the output file
    }

    // Process the miRNA file, filter out rows where the miRNA targets an mRNA that is also targeted by RBP
    int i = 0;
    while (getline(in_mirna, line)) {
        vector<string> fields = parse_line(line);
        if (fields.size() < 7) continue;  // Skip incomplete rows

        string s7 = fields[6];  // miRNA target gene

        // Write to output file if the target gene is not already in the RBP set
        if (rbp_genes.find(s7) == rbp_genes.end()) {
            out_mirna << line << endl;
        }

        i++;
        if (i % 100 == 0) cout << "Processed " << i << " miRNA entries." << endl;
    }

    cout << "Merging and filtering completed." << endl;
}

// Merge RBP and miRNA interaction data based on mRNA
void merge_rbp_mirna_interactions(const string& rbp_file, const string& mirna_file, const string& output_file) {
    map<string, vector<string>> rbp_map;  // To store RBP-mRNA interaction data
    ifstream in_rbp(rbp_file);
    ifstream in_mirna(mirna_file);
    ofstream out(output_file);

    string line;

    // Check if the input files were successfully opened
    if (!in_rbp || !in_mirna) {
        cerr << "Error: Unable to open input files." << endl;
        return;
    }

    // Read the RBP file and store the data in a map
    while (getline(in_rbp, line)) {
        stringstream ss(line);
        string rbp, mRNA, extra_info, start, end, strand, gene;
        ss >> rbp >> mRNA >> extra_info >> start >> end >> strand >> gene;
        string combined_info = start + '\t' + gene + '\t' + rbp + '\t' + mRNA + '\t' + extra_info;
        rbp_map[gene].push_back(combined_info);  // Store the data by gene
    }

    cout << "Step 1: Merging miRNA-mRNA and RBP-mRNA interactions..." << endl;

    int line_count = 0;
    // Read the miRNA file and check if there are corresponding RBP interactions
    while (getline(in_mirna, line)) {
        stringstream ss(line);
        string mirna, mRNA, extra_info, start, end, strand, gene;
        ss >> mirna >> mRNA >> extra_info >> start >> end >> strand >> gene;
        string combined_info = start + '\t' + gene + '\t' + mirna + '\t' + mRNA + '\t' + extra_info;
        
        // If the miRNA target mRNA exists in the RBP interactions, write it to the output
        if (rbp_map.count(gene)) {
            for (const auto& rbp_info : rbp_map[gene]) {
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

// Analyze overlapping and adjacent interactions between RBP and miRNA on mRNA
void analyze_interactions(const string& input_file, const string& output_file) {
    ifstream in(input_file);
    ofstream out(output_file);
    string line;

    // Maps to store total, overlap, and adjacent interaction counts
    map<string, int> total_count, overlap_count, adjacent_count;

    if (!in) {
        cerr << "Error: Unable to open the input file." << endl;
        return;
    }

    int line_count = 0;
    while (getline(in, line)) {
        stringstream ss(line);
        string rbp, mirna, mRNA_rbp, mRNA_mirna, gene_rbp, gene_mirna;
        int start_rbp, end_rbp, start_mirna, end_mirna;
        
        ss >> rbp >> gene_rbp >> mRNA_rbp >> start_rbp >> end_rbp 
           >> mirna >> gene_mirna >> mRNA_mirna >> start_mirna >> end_mirna;

        string interaction_key = rbp + '\t' + gene_rbp + '\t' + mirna;  // Unique RBP-miRNA pair

        total_count[interaction_key]++;  // Count total interactions

        // Check if interactions are overlapping or adjacent
        if (mRNA_rbp == mRNA_mirna && 
            abs(end_rbp - start_rbp) + abs(end_mirna - start_mirna) > 
            max({abs(end_rbp - start_mirna), abs(end_rbp - end_mirna), abs(start_rbp - start_mirna), abs(start_rbp - end_mirna)})) {
            overlap_count[interaction_key]++;
        } else {
            adjacent_count[interaction_key]++;
        }

        line_count++;
        if (line_count % 100000 == 0) {
            cout << "Processed " << line_count << " interactions." << endl;
        }
    }

    // Output the analysis results
    for (const auto& entry : total_count) {
        out << entry.first << '\t' << entry.second << '\t' 
            << overlap_count[entry.first] << '\t' << adjacent_count[entry.first] << endl;
    }

    cout << "Analysis completed." << endl;
}

int main() {
    string rbp_file = "./inputs/RBPmRNAv27.bed";
    string mirna_file = "./inputs/miRNAmRNAv27.bed";
    string intermediate_file1 = "./outputs/RBPmRNAv27_filtered.bed";
    string intermediate_file2 = "./outputs/miRNAmRNAv27_filtered.bed";
    string merged_file = "./outputs/RBP_m_mi_intermediate.bed";
    string output_file = "./outputs/RBP_m_mi_analysis.bed";

    // Step 1: Filter RBP and miRNA based on mRNA gene
    cout << "Step 1: Merging and filtering RBP and miRNA based on mRNA gene..." << endl;
    merge_mrna_interactions(rbp_file, mirna_file, intermediate_file1, intermediate_file2);
    
    // Step 2: Merge RBP and miRNA interaction data
    cout << "Step 2: Merging filtered RBP and miRNA interactions..." << endl;
    merge_rbp_mirna_interactions(intermediate_file1, intermediate_file2, merged_file);

    // Step 3: Analyze adjacent and overlapping interactions
    cout << "Step 3: Analyzing adjacent and overlapping interactions..." << endl;
    analyze_interactions(merged_file, output_file);

    return 0;
}
