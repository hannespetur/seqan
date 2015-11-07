#ifndef SEQAN_HTS_IO_HTS_FILE_H_
#define SEQAN_HTS_IO_HTS_FILE_H_

#include <seqan/hts_io/hts_alignment_record.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>


namespace seqan {

class HtsFile
{
  public:
    const char * filename;
    htsFile * file_pointer;
    bam_hdr_t * hdr;
    bam1_t * hts_record;

    HtsFile(const char * filename)
        : filename(filename)
        {}
};

inline bool
readRecord(HtsSequenceRecord & record, HtsFile & file)
{
    if (sam_read1(file.file_pointer, file.hdr, file.hts_record) < 0)
        return false;

    record.qName = bam_get_qname(file.hts_record);
    int32_t lqseq = file.hts_record->core.l_qseq;
    resize(record.seq, lqseq);
    uint8_t* seqptr = bam_get_seq(file.hts_record);

    for (int i = 0; i < lqseq; ++i)
    {
      record.seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }
    return true;
}

inline void
open(HtsFile & file)
{
    std::cout << "filename = " << file.filename << std::endl;
    file.file_pointer = hts_open(file.filename, "r");
    if (file.file_pointer == NULL)
    {
        SEQAN_FAIL("Could not open file with filename %s", file.filename);
    }
    
    file.hdr = sam_hdr_read(file.file_pointer);
    file.hts_record = bam_init1();
}

// inline void
// loadIndex(HtsFile & file)
// {
//     hts_idx_t * sam_index_load(file, const char *fn);
// }

} // namespace seqan


#endif  // SEQAN_HTS_IO_HTS_FILE_H_
