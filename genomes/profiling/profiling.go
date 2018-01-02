package profiling

import (
	"fmt"
	"os"

	"github.com/kussell-lab/biogo/feat/gff"
	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/ncbiftp/seqrecord"
	"github.com/kussell-lab/ncbiftp/taxonomy"
)

// Pos contains genomic position profile.
type Pos struct {
	Base byte
	Type byte
	Gene string
}

// A position could be one of those:
// 	First site in a codon;
// 	Second site in a codon;
//	Third site in a codon;
// 	Fourfold site
// 	NonCoding site
// 	Undefined
const (
	NonCoding byte = '0'
	FirstPos  byte = '1'
	SecondPos byte = '2'
	ThirdPos  byte = '3'
	FourFold  byte = '4'
	Undefined byte = '5'
	Coding    byte = '6'
)

func ProfileGenome(genome []byte, gffRecords []*gff.Record, gc *taxonomy.GeneticCode) (profile []Pos) {

	// mark all sites as non-coding.
	profile = make([]Pos, len(genome))
	for i := 0; i < len(profile); i++ {
		profile[i] = Pos{Type: NonCoding}
	}

	// for each gene, mark codon positions.
	geneIndex := 0
	for _, rec := range gffRecords {
		geneIndex++
		// prepare nucleotide sequence,
		// we need it for determine 4-fold codons.
		var nucl []byte
		if rec.End >= rec.Start {
			nucl = genome[rec.Start-1 : rec.End]
		} else {
			// skip genes across boundary.
			continue
		}

		// reverse and complement the negative strain.
		if rec.Strand == gff.ReverseStrand {
			nucl = seq.Complement(seq.Reverse(nucl))
		}

		prof := make([]byte, len(nucl))
		for j, _ := range nucl {
			switch (j + 1) % 3 {
			case 1:
				prof[j] = FirstPos
			case 2:
				prof[j] = SecondPos
			case 0:
				// determine if it is a fourfold site.
				codon := nucl[j-2 : j+1]
				if gc.FFCodons[string(codon)] {
					prof[j] = FourFold
				} else {
					prof[j] = ThirdPos
				}
			}
		}

		// if it is a negative strain, reverse the profile to match the positive strain.
		if rec.Strand == gff.ReverseStrand {
			prof = seq.Reverse(prof)
		}

		// write the position profile into the entire genomic profile.
		for j, p := range prof {
			index := rec.Start - 1 + j
			// check overlapping.
			// if overlap, simply mark it as undefined.
			base := genome[index]
			if profile[index].Type == NonCoding {
				profile[index] = Pos{Type: p, Base: base, Gene: fmt.Sprintf("%s_%d", rec.SeqName, geneIndex)}
			} else {
				profile[index] = Pos{Type: Undefined, Base: base, Gene: fmt.Sprintf("%s_%d", rec.SeqName, geneIndex)}
			}
		}
	}

	return
}

// Generate codon position profile for the entire genome.
// First, we mark every position as NonCoding.
// Then, for each coding (gene) region, we determine each codon position.
// If there is an overlapping region between two genes, mark them as undefined.
func ProfileGenome1(genomeFileName, pttFileName string, gc *taxonomy.GeneticCode) (profile []Pos) {
	// read .ptt file and obtain gene coding region.
	ptts := readPtt(pttFileName)
	// read genome sequence.
	genome := readGenome(genomeFileName)
	s := genome.Seq

	// mark all sites as non-coding.
	profile = make([]Pos, len(s))
	for i := 0; i < len(profile); i++ {
		profile[i] = Pos{Type: NonCoding}
	}

	// for each gene, mark codon positions.
	for _, ptt := range ptts {
		// prepare nucleotide sequence,
		// we need it for determine 4-fold codons.
		var nucl []byte
		if ptt.Loc.To >= ptt.Loc.From {
			nucl = s[ptt.Loc.From-1 : ptt.Loc.To]
		} else {
			// skip genes across boundary.
			continue
		}

		// reverse and complement the negative strain.
		if ptt.Loc.Strand == "-" {
			nucl = seq.Complement(seq.Reverse(nucl))
		}

		prof := make([]byte, len(nucl))
		for j, _ := range nucl {
			switch (j + 1) % 3 {
			case 1:
				prof[j] = FirstPos
			case 2:
				prof[j] = SecondPos
			case 0:
				// determine if it is a fourfold site.
				codon := nucl[j-2 : j+1]
				if gc.FFCodons[string(codon)] {
					prof[j] = FourFold
				} else {
					prof[j] = ThirdPos
				}
			}
		}

		// if it is a negative strain, reverse the profile to match the positive strain.
		if ptt.Loc.Strand == "-" {
			prof = seq.Reverse(prof)
		}

		// write the position profile into the entire genomic profile.
		for j, p := range prof {
			index := ptt.Loc.From - 1 + j
			// check overlapping.
			// if overlap, simply mark it as undefined.
			base := s[index]
			if profile[index].Type == NonCoding {
				profile[index] = Pos{Type: p, Base: base, Gene: ptt.PID}
			} else {
				profile[index] = Pos{Type: Undefined, Base: base, Gene: ptt.PID}
			}
		}
	}

	return
}

// read ptt file.
func readPtt(fileName string) []seqrecord.Ptt {
	reader := seqrecord.NewPttFile(fileName)
	ptts := reader.ReadAll()
	return ptts
}

// read genome sequence from a FASTA file.
func readGenome(fileName string) *seq.Sequence {
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	reader := seq.NewFastaReader(f)
	sequences, err := reader.ReadAll()
	if err != nil {
		panic(err)
	}

	return sequences[0]
}
