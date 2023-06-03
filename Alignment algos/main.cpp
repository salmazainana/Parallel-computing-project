#include <iostream>
#include <vector>
#include <algorithm>
#include <thread>
#include <chrono>

#include "Mmalgorithm.cpp"
#include "gotoh.cpp"

#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
#include <curl/curl.h>
#include <rapidjson/document.h>

using namespace std;

// The implementation of the UniProtKB API is given below

// Struct to hold the data received from the HTTP response
struct HttpResponse {
    string data;
    int code;
};

// Callback function to write HTTP response data
size_t WriteCallback(void* contents, size_t size, size_t nmemb, string* output) {
    size_t totalSize = size * nmemb;
    output->append((char*)contents, totalSize);
    return totalSize;
}

// Function to perform an HTTP GET request
HttpResponse performHttpRequest(const string& url) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        cerr << "Failed to initialize CURL." << endl;
        return { "", -1 };
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    string responseBuffer;
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &responseBuffer);

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        cerr << "Failed to perform HTTP request: " << curl_easy_strerror(res) << endl;
        return { "", -1 };
    }

    long httpCode = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &httpCode);

    curl_easy_cleanup(curl);

    return { responseBuffer, static_cast<int>(httpCode) };
}

// Function to retrieve sequence from UniProtKB based on accession
string getSequenceFromUniProtKB(const string& accession) {
    string url = "https://www.uniprot.org/uniprot/" + accession + ".fasta";
    HttpResponse response = performHttpRequest(url);
    if (response.code != 200) {
        cerr << "Failed to retrieve sequence from UniProtKB for accession: " << accession << endl;
        return "";
    }

    // Parse the FASTA response to extract the sequence
    size_t sequenceStart = response.data.find('\n');
    if (sequenceStart != string::npos && sequenceStart < response.data.length() - 1) {
        return response.data.substr(sequenceStart + 1);
    }

    cerr << "Failed to parse sequence from UniProtKB response for accession: " << accession << endl;
    return "";
}

int main() {
    string accession1 = "P31964";
    string accession2 = "Q62258";

    string seq1 = getSequenceFromUniProtKB(accession1);
    string seq2 = getSequenceFromUniProtKB(accession2);

    if (seq1.empty() || seq2.empty()) {
        cerr << "Failed to retrieve sequences from UniProtKB." << endl;
        return 1;
    }

    auto begin = chrono::high_resolution_clock::now();
    sequenceAlignmentParallel(seq1, seq2, 8);
    auto end = chrono::high_resolution_clock::now();
    cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;

    auto begin = chrono::high_resolution_clock::now();
    gotohAlignmentParallel(seq1, seq2, 8);
    auto end = chrono::high_resolution_clock::now();
    cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms by Gotoh Algorithm" << endl;

    //auto begin = chrono::high_resolution_clock::now();
    //NWParallel(seq1, seq2, 8);
    //auto end = chrono::high_resolution_clock::now();
    //cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms by Needlman-Wunch Algorithm" << endl;

    return 0;
}
