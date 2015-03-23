package reports

import (
	"bufio"
	"io"
	"strings"
)

// To retreive prokaryotes information
// from the folder of GENOME_REPORTS in
// ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS

// Strain information.
type Strain struct {
	ProjectId   string   // BioProject ID.
	Name        string   // Strain name.
	TaxId       string   // Taxonomy ID.
	Genomes     []Genome // Genome RefSeq accessions.
	Path        string   // folder path to the NCBI ftp.
	GeneticCode string   // Genetic codon table ID.
	Status      string   // Status, complete or not.
}

type Genome struct {
	Accession  string
	Replicon   string
	Length     int
	Seq        []byte
	PosProfile []byte
}

// Read prokaryotes.txt
func ReadProkaryotes(f io.Reader) (strains []Strain) {
	// create a buffer reader.
	rd := bufio.NewReader(f)

	// read the first commented line to
	// determine the field names.
	nameMap := make(map[string]int)
	if r1, _, err := rd.ReadRune(); err == nil {
		if r1 == '#' {
			line, err := rd.ReadString('\n')
			if err != nil {
				panic(err)
			} else {
				names := strings.Split(strings.TrimSpace(line), "\t")
				for i := 0; i < len(names); i++ {
					nameMap[names[i]] = i
				}
			}
		}
	} else {
		panic(err)
	}

	records := [][]string{}
	for {
		line, err := rd.ReadString('\n')
		// continue, if it is a comment line.
		if line[0] == '#' {
			continue
		}

		if err != nil {
			if err != io.EOF {
				panic(err)
			} else {
				break
			}
		} else {
			fields := strings.Split(strings.TrimSpace(line), "\t")
			records = append(records, fields)
		}
	}

	for _, fields := range records {
		s := Strain{}
		s.Name = fields[nameMap["Organism/Name"]]
		s.TaxId = fields[nameMap["TaxID"]]
		s.ProjectId = fields[nameMap["BioProject ID"]]
		s.Path = fields[nameMap["FTP Path"]]
		s.Status = fields[nameMap["Status"]]

		chromosomes := fields[nameMap["Chromosomes/RefSeq"]]
		// remove redundant.
		m := make(map[string]bool)
		for _, g := range strings.Split(chromosomes, ",") {
			acc := strings.Split(strings.TrimSpace(g), ".")[0]
			m[acc] = true
		}

		for acc, _ := range m {
			s.Genomes = append(s.Genomes,
				Genome{Accession: acc, Replicon: "Chromosome"})
		}

		strains = append(strains, s)
	}

	return
}
