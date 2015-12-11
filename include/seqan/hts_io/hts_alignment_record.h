#ifndef SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_
#define SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_

#include <seqan/basic.h>

#include <cstdio>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

namespace seqan {

class HtsSequenceRecord
{
  public:
    String<char> qName;
    String<Iupac> seq;

    static __int32 const INVALID_POS = -1;
    static __int32 const INVALID_REF_ID = -1;
    static __int32 const INVALID_LEN = 0;
    static __uint32 const INVALID_QID = 4294967295u;

    HtsSequenceRecord()
      : qName(), seq() {}

    HtsSequenceRecord(bam1_t * hts_record)
    {
        HtsSequenceRecord::parse(hts_record);
    }

    virtual void parse(bam1_t * hts_record)
    {
        qName = bam_get_qname(hts_record);
        int32_t lqseq = hts_record->core.l_qseq;
        resize(seq, lqseq, Exact());
        uint8_t* seqptr = bam_get_seq(hts_record);

        for (int i = 0; i < lqseq; ++i)
        {
          seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
        }
    }
};


class HtsAlignmentRecord : public HtsSequenceRecord
{
  public:
    // __uint32 _qID; What is this?
    // uint16_t flag;
    StringSet<String<char> > cigar; // TODO: Change to String<CigarElement>?
    String<char> qual;
    // String<char> tags;

    void parse(bam1_t * hts_record)
    {
        HtsSequenceRecord::parse(hts_record);

        // parse quality
        uint8_t* qualptr = bam_get_qual(hts_record);
        resize(qual, length(seq), Exact());

        for (unsigned i = 0; i < length(seq); ++i, ++qualptr)
        {
            qual[i] = static_cast<char>(*qualptr + 33);
        }

        // parse cigar
        uint32_t* cigarptr = bam_get_cigar(hts_record);

        for (unsigned i = 0; i < (hts_record)->core.n_cigar; ++i, ++cigarptr)
        {
            static char const * CIGAR_MAPPING = "MIDNSHP=X*******";
            char buffer [16];
            sprintf (buffer, "%c%d ", static_cast<char>(CIGAR_MAPPING[(*cigarptr >> 28) & 0xF]), (*cigarptr >> 4) & 0xFFFFF);
            appendValue(cigar, buffer);
        }
    }
};

} // namespace seqan

#endif // SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_
