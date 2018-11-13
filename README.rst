SeqAn - The Library for Sequence Analysis (now with HTS!)
=========================================================

This is a fork of SeqAn 2.1 that uses htslib to read and write VCF/BCF and SAM/BAM/CRAM files. List of new added features:
 * VCF/BCF file support.
    - Reading implemented.
    - Writing **not** implemented.
    - Indexes (Tabix)
        + Reading implemented.
        + Building **not** implemented.
    - Parsing to SeqAn's VCF record is implemented.
 * HTS (SAM/BAM/CRAM) support.
    - Reading implemented.
    - Writing implemented.
    - Indexes (.bai, .csi, .crai)
        + Reading implemented.
        + Building implemented.
    - Parsing to SeqAn alignment records.


If you need a feature which is not implemented you can ask me (Hannes Pétur, hannese@decode.is) to add it, or use the `htslib API <https://github.com/samtools/htslib>`_ directly.


Compilation differences between SeqAn and SeqAnHTS
--------------------------------------------------

Use the required SeqAn compilation flags and also include htslib by setting ``-I<path_to_htslib/include>`` and link to the htslib library.


API differences between SeqAn and SeqAnHTS
------------------------------------------

The 'BamIndex' object
~~~~~~~~~~~~~~~~~~~~~

BamIndex and BamFileIn (or BamFileOut) have been merged into one class called BamFileIn (or BamFileOut). A pointer to the index is simply an instance variable of BamFile. If no index has been loaded or created it is a null pointer. For example, this:

.. code-block:: cpp

  BamFileIn bamFileIn;
  open(bamFileIn, ...);
  CharString baiPathIn = "my/path/to/index";
  BamIndex<Bai> baiIndex;
  if (!open(baiIndex, toCString(baiPathIn)))
  {
    std::cerr << "ERROR: Could not read BAI index file " << baiPathIn << "\n";
    return 1;
  }

can be changed to this

.. code-block:: cpp

  BamFileIn bamFileIn;
  open(bamFileIn, ...);
  CharString baiPathIn = "my/path/to/index";
  if (!loadIndex(bamFileIn, toCString(baiPathIn)))
  {
    std::cerr << "ERROR: Could not read BAI index file " << baiPathIn << "\n";
    return 1;
  }


The 'context' function
~~~~~~~~~~~~~~~~~~~~~~
The 'context' function is not defined or any other function that uses it for two reasons: One, there is nothing in htslib which is a sensible replacement object to it, and two, this is probably something the user doesn't need to worry about. A typical example is if the user wants to read a specified region of an alignment file. In that case the user must first use the context to figure out which rID the was specified for the region in the header.

.. code-block:: cpp

  CharString chr = "chr19";
  int start = 3000;
  int end = 4000;
  int rID = 0;
  if (!getIdByName(rID, contigNamesCache(context(bamFileIn)), chr))
  {
    std::cerr << "ERROR: Reference sequence named " << chr << " not known.\n";
    return 1;
  }
  bool hasAlignments = false;
  if (!jumpToRegion(bamFileIn, hasAlignments, rID, start, end, baiIndex))
  {
    std::cerr << "ERROR: Could not jump to " << start << ":" << end << "\n";
    return 1;
  }

This can be changed to

.. code-block:: cpp

  CharString chr = "chr19";
  int start = 3000;
  int end = 4000;
  if (!setRegion(bamFileIn, toCString(chr), start, end))
  {
    std::cerr << "ERROR: Could not jump to " << chr << ":" << start << "-" << end << "\n";
    return 1;
  }


Examples
--------
BCF/Tabix example
~~~~~~~~~~~~~~~~~
.. code-block:: cpp

  seqan::Tabix index;
  seqan::open(index, "/path/to/my/file/example.vcf.gz");
  seqan::setRegion(index, "chrX:A-B"); // "chrX" and "chrX:A" also supported
  seqan::VcfRecord record;

  while (seqan::readRegion(record, index))
  {
    // Do stuff with record
  }

HTS file read and write example (SAM/BAM/CRAM format is automatically detected)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

  seqan::HtsFile read_file("/path/to/some/existing/file.cram", "r");
  seqan::HtsFile write_file("/path/to/a/new/file.cram", "wb"); // binary mode required for BAM and CRAM

  seqan::copyHeader(write_file, read_file);
  seqan::writeHeader(write_file);

  while (seqan::readRecord(read_file))
  {
    seqan::copyRecord(write_file, read_file);
    // Here you could change the write_file.hts_record if you want.
    seqan::writeRecord(write_file);
  }
  // Here we have copied all records of read_file and written them to write_file


HTS file index example
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: cpp

  seqan::HtsFile hts_file("/path/to/some/existing/file.cram", "r");

  if (!seqan::loadIndex(hts_file))
  {
    // Build it if we cannot find it
    seqan::buildIndex(hts_file);
    seqan::loadIndex(hts_file);
  }

  seqan::setRegion(hts_file, "chrX:A-B");
  seqan::HtsSequenceRecord record; // Only parses qName and sequence, use seqan::BamAlignmentRecord to parse all

  while (readRegion(record, hts_file))
  {
    // Do stuff with each record that overlaps the chrX:A-B region.
  }
