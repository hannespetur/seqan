#ifndef SEQAN_BAM_IO_CRAM_FILE_H_
#define SEQAN_BAM_IO_CRAM_FILE_H_

#include <seqan/bam_io/cram_alignment_record.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>


namespace seqan {

class CramFile
{
  public:
    samFile *in;
    bam_hdr_t *hdr;
    bam1_t *b;

    ~CramFile()
    {
        if (in == NULL)
            return;

        bam_destroy1(b);
        bam_hdr_destroy(hdr);
        hts_close(in);
    }
};

inline bool
readRecord(CramSequenceRecord & record, CramFile & file)
{
    record.qName = bam_get_qname(file.b);

    int32_t lqseq = file.b->core.l_qseq;
    resize(record.seq, lqseq);
    uint8_t* seqptr = bam_get_seq(file.b);

    for (int i = 0; i < lqseq; ++i)
    {
      record.seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }
}

inline void
open(CramFile & file, const char * filename)
{
    file.in = hts_open(filename, "r");
    if (file.in == NULL) return;
    file.hdr = sam_hdr_read(file.in);
    file.b = bam_init1();
}

} // namespace seqan


#endif  // SEQAN_BAM_IO_CRAM_FILE_H_
