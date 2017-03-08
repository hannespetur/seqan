#include <seqan/basic.h>
#include <seqan/vcf_io.h>

int main(int argc, char* argv[]) {
    seqan::Tabix index;
    seqan::CharString h_string;
    seqan::open(index, "/odinn/data/results/seqgts/aster/vcf/chr3/005549890-005599888.raw.vcf.gz");
    seqan::setRegion(index, "chr3:5549990");
    seqan::getHeader(h_string, index);
    seqan::VcfRecord record;

    while (seqan::readRegion(record, index)) {
        std::cout << seqan::length(index.samples) << " - " << seqan::length(record.genotypeInfos) << std::endl;
        std::cout << index.samples[0] << ": " << record.genotypeInfos[0] << std::endl;
    }
}

