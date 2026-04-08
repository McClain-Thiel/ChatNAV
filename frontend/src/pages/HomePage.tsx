import { useNavigate } from 'react-router-dom';
import SectionRule from '../components/SectionRule';

export default function HomePage() {
  const navigate = useNavigate();

  return (
    <div style={{
      height: '100%',
      overflowY: 'auto',
      background: 'var(--bg-page)',
    }}>
      <div style={{ maxWidth: 720, margin: '0 auto', padding: '48px 40px 80px' }}>
        {/* Hero */}
        <div style={{ marginBottom: 40 }}>
          <div style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 9,
            letterSpacing: 2.5,
            textTransform: 'uppercase',
            color: 'var(--text-muted)',
            marginBottom: 8,
          }}>
            Computational Vaccine Design
          </div>
          <h1 style={{
            fontFamily: "'EB Garamond', serif",
            fontSize: 42,
            fontWeight: 500,
            color: 'var(--text-primary)',
            lineHeight: 1.15,
            marginBottom: 12,
            letterSpacing: -0.5,
          }}>
            Personalised mRNA Neoantigen Vaccines
          </h1>
          <p style={{
            fontFamily: "'EB Garamond', serif",
            fontStyle: 'italic',
            fontSize: 19,
            color: 'var(--text-meta)',
            lineHeight: 1.5,
            marginBottom: 24,
          }}>
            From somatic mutations to a synthesis-ready mRNA sequence in a single pipeline.
          </p>
          <div style={{ borderTop: '2px solid var(--border-heavy)' }} />
          <div style={{ borderTop: '1px solid var(--border-mid)', marginTop: 4 }} />
        </div>

        {/* What it does */}
        <div style={{ marginBottom: 36 }}>
          <p style={bodyStyle}>
            ChatNAV takes a patient's somatic mutation data (MAF or VCF), HLA types,
            and gene expression profile as input and produces a complete, synthesis-ready
            mRNA vaccine sequence as output. The pipeline identifies the mutations most
            likely to trigger an immune response, assembles them into an optimised
            polyepitope construct, and designs the mRNA with codon and structural
            optimisation for maximum translational efficiency.
          </p>
        </div>

        {/* Pipeline overview */}
        <SectionRule label="How It Works" />
        <div style={{ marginBottom: 36 }}>
          {PIPELINE_STEPS.map((step, i) => (
            <div key={i} style={{
              display: 'flex',
              gap: 16,
              padding: '12px 0',
              borderBottom: '1px solid var(--border-faint)',
            }}>
              <div style={{
                fontFamily: "'Inconsolata', monospace",
                fontSize: 9,
                color: 'var(--text-faint)',
                minWidth: 20,
                marginTop: 2,
              }}>
                {String(i + 1).padStart(2, '0')}
              </div>
              <div style={{ flex: 1 }}>
                <div style={{
                  fontFamily: "'EB Garamond', serif",
                  fontSize: 15,
                  fontWeight: 600,
                  color: 'var(--text-primary)',
                  marginBottom: 2,
                }}>
                  {step.title}
                </div>
                <div style={{
                  fontFamily: "'Spectral', serif",
                  fontSize: 13.5,
                  lineHeight: 1.7,
                  color: 'var(--text-secondary)',
                }}>
                  {step.description}
                </div>
              </div>
              <div style={{
                fontFamily: "'Inconsolata', monospace",
                fontSize: 8.5,
                color: 'var(--text-muted)',
                whiteSpace: 'nowrap',
                marginTop: 2,
              }}>
                {step.tool}
              </div>
            </div>
          ))}
        </div>

        {/* Example output */}
        <SectionRule label="Example Output" />
        <div style={{
          background: 'var(--bg-inset)',
          border: '1px solid var(--border-light)',
          borderRadius: 3,
          padding: '16px 20px',
          marginBottom: 36,
          fontFamily: "'Inconsolata', monospace",
          fontSize: 10,
          lineHeight: 1.8,
          color: 'var(--text-secondary)',
        }}>
          <div>Patient: TCGA-EE-A3J5 (TCGA SKCM melanoma)</div>
          <div>Variants: 2,584 somatic → 1,890 peptide windows → 260 strong binders</div>
          <div>&nbsp;&nbsp;→ 84 pass filters → 20 selected</div>
          <div>Polyepitope: 322 amino acids</div>
          <div>mRNA: 1,253 nt | GC: 57.2% | MFE: −679.9 kcal/mol | CAI: 0.900</div>
        </div>

        <button
          onClick={() => navigate('/jobs/demo')}
          style={{
            display: 'block',
            width: '100%',
            padding: '14px 0',
            fontFamily: "'EB Garamond', serif",
            fontSize: 15,
            fontWeight: 500,
            color: 'var(--text-primary)',
            background: 'transparent',
            border: '1px solid var(--border-heavy)',
            borderRadius: 2,
            cursor: 'pointer',
            marginBottom: 40,
          }}
        >
          View Demo Results →
        </button>

        {/* Important caveats */}
        <SectionRule label="Important Caveats" />
        <div style={{
          borderLeft: '2px solid var(--status-running)',
          paddingLeft: 18,
          marginBottom: 36,
        }}>
          <p style={{ ...bodyStyle, marginBottom: 16 }}>
            Neoantigen vaccines are an area of active clinical research with
            promising but preliminary evidence. Users should understand the following
            limitations before interpreting pipeline outputs:
          </p>
          <ul style={{ ...bodyStyle, paddingLeft: 18, marginBottom: 0 }}>
            <li style={{ marginBottom: 10 }}>
              <strong>Not a clinical tool.</strong> This pipeline is for research purposes
              only. Outputs have not been validated in a clinical setting and should not
              be used to make treatment decisions without appropriate regulatory review.
            </li>
            <li style={{ marginBottom: 10 }}>
              <strong>Prediction uncertainty.</strong> Immunogenicity prediction remains
              an unsolved problem. The best available models (BigMHC-IM) achieve AUC ~0.92
              on benchmark data, but real-world performance varies by tumour type, HLA
              context, and microenvironment.
            </li>
            <li style={{ marginBottom: 10 }}>
              <strong>The probability estimate is a model output, not a guarantee.</strong> The
              "P(≥1 active)" metric uses a simplified binomial model: P = 1 − (1 − p)<sup>n</sup>.
              It assumes epitopes are independent and does not account for immune escape,
              tolerance, or tumour heterogeneity.
            </li>
            <li style={{ marginBottom: 10 }}>
              <strong>CCF limitations.</strong> Cancer cell fraction estimation from
              MAF-only input defaults to 1.0 (all mutations assumed clonal). This
              overstates the clonality of subclonal mutations. BAM + CNV data enables
              accurate CCF estimation via PyClone-VI.
            </li>
            <li style={{ marginBottom: 10 }}>
              <strong>Expression approximation.</strong> Without patient RNA-seq,
              expression is imputed from GTEx tissue-type medians — a population average,
              not the patient's actual tumour expression.
            </li>
            <li style={{ marginBottom: 10 }}>
              <strong>HLA Class II gaps.</strong> MHC Class II prediction (MHCnuggets) is
              less accurate than Class I (MHCflurry). CD4+ T cell epitopes are important
              for sustained immune responses but harder to predict.
            </li>
          </ul>
        </div>

        {/* Scoring methodology */}
        <SectionRule label="Scoring Methodology" />
        <div style={{ marginBottom: 36 }}>
          <p style={{ ...bodyStyle, marginBottom: 16 }}>
            Each candidate neoantigen receives a composite score combining multiple
            orthogonal signals, weighted to balance immunogenicity prediction with
            practical vaccine design constraints:
          </p>
          <div style={{
            display: 'grid',
            gridTemplateColumns: '1fr 1fr',
            gap: '8px 24px',
            marginBottom: 16,
          }}>
            {SCORING_COMPONENTS.map(comp => (
              <div key={comp.name} style={{
                padding: '8px 0',
                borderBottom: '1px solid var(--border-faint)',
              }}>
                <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                  <span style={{
                    fontFamily: "'Inconsolata', monospace",
                    fontSize: 10,
                    color: 'var(--text-secondary)',
                  }}>
                    {comp.name}
                  </span>
                  <span style={{
                    fontFamily: "'Inconsolata', monospace",
                    fontSize: 10,
                    color: 'var(--text-muted)',
                  }}>
                    {comp.weight}
                  </span>
                </div>
                <div style={{
                  fontFamily: "'Spectral', serif",
                  fontSize: 11.5,
                  color: 'var(--text-meta)',
                  marginTop: 2,
                }}>
                  {comp.desc}
                </div>
              </div>
            ))}
          </div>

          <p style={{
            fontFamily: "'Spectral', serif",
            fontStyle: 'italic',
            fontSize: 12,
            color: 'var(--text-faint)',
            lineHeight: 1.6,
          }}>
            Hard filters: binding rank ≤ 2%, TPM ≥ 1.0, CCF ≥ 0.5.
            Bonuses: frameshift (+0.10), shared neoantigen (+0.15), CD4 epitope (+0.05).
          </p>
        </div>

        {/* References */}
        <SectionRule label="Key References" />
        <div style={{ marginBottom: 20 }}>
          {REFERENCES.map((ref, i) => (
            <div key={i} style={{
              fontFamily: "'Spectral', serif",
              fontStyle: 'italic',
              fontSize: 12,
              color: 'var(--text-meta)',
              marginBottom: 6,
              lineHeight: 1.5,
            }}>
              <span style={{ fontStyle: 'normal', color: 'var(--text-muted)', marginRight: 6 }}>
                [{i + 1}]
              </span>
              {ref}
            </div>
          ))}
        </div>

        {/* Footer */}
        <div style={{ borderTop: '2px solid var(--border-heavy)', paddingTop: 16 }}>
          <div style={{
            fontFamily: "'Spectral', serif",
            fontStyle: 'italic',
            fontSize: 12,
            color: 'var(--text-faint)',
            textAlign: 'center',
          }}>
            ChatNAV v0.1.0 · Pipeline code: MIT · External tools have their own licences
          </div>
        </div>
      </div>
    </div>
  );
}

const bodyStyle: React.CSSProperties = {
  fontFamily: "'Spectral', serif",
  fontSize: 15,
  lineHeight: 1.85,
  color: 'var(--text-secondary)',
  textAlign: 'justify',
};

const PIPELINE_STEPS = [
  { title: 'Alignment', description: 'Map paired-end reads to GRCh38 reference genome.', tool: 'BWA-MEM2' },
  { title: 'Variant Calling', description: 'Identify somatic SNVs and indels; annotate with VEP.', tool: 'Mutect2 + VEP' },
  { title: 'HLA Typing', description: 'Determine patient HLA class I and II alleles from normal DNA.', tool: 'HLA-HD' },
  { title: 'Expression', description: 'Quantify gene expression from RNA-seq or GTEx tissue priors.', tool: 'Salmon' },
  { title: 'Clonality', description: 'Estimate cancer cell fraction per mutation; prioritise truncal variants.', tool: 'PyClone-VI' },
  { title: 'Candidate Generation', description: 'Generate 8-11mer peptide windows around each somatic mutation.', tool: 'pVACseq' },
  { title: 'MHC Binding', description: 'Predict peptide-HLA binding affinity for the patient\'s alleles.', tool: 'MHCflurry 2.0' },
  { title: 'Immunogenicity', description: 'Score with BigMHC-IM immunogenicity, proteome foreignness, and agretopicity.', tool: 'BigMHC-IM' },
  { title: 'Ranking', description: 'Composite scoring with configurable weights; select top 20 with HLA diversity.', tool: 'Weighted composite' },
  { title: 'Polyepitope', description: 'Assemble epitopes with optimised ordering, linkers, and signal peptide.', tool: 'Greedy TSP' },
  { title: 'mRNA Design', description: 'Codon and structure optimisation; add UTRs, poly-A, pseudouridine positions.', tool: 'LinearDesign' },
];

const SCORING_COMPONENTS = [
  { name: 'Immunogenicity', weight: '0.35', desc: 'BigMHC-IM transformer prediction' },
  { name: 'Foreignness', weight: '0.15', desc: 'k-mer distance from human proteome' },
  { name: 'Agretopicity', weight: '0.10', desc: 'log2(wt_rank / mut_rank)' },
  { name: 'Binding', weight: '0.15', desc: 'MHCflurry percentile rank' },
  { name: 'Stability', weight: '0.10', desc: 'pMHC complex half-life' },
  { name: 'Expression', weight: '0.10', desc: 'Normalised TPM from RNA-seq' },
  { name: 'Clonality (CCF)', weight: '0.05', desc: 'Cancer cell fraction' },
  { name: 'Structural', weight: '0.00–0.05', desc: 'TCR exposure (tiered: position/PANDORA/AlphaFold)' },
];

const REFERENCES = [
  'Albert et al. (2023). BigMHC: accurate prediction of MHC-I antigen presentation and immunogenicity. Nature Machine Intelligence.',
  'O\'Donnell et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides. Cell Systems.',
  'Zhang et al. (2023). LinearDesign: Efficient algorithms for optimized mRNA sequence design. Nature.',
  'Gartner et al. (2021). A machine learning model for ranking candidate HLA-I neoantigens. Nature Medicine.',
  'Hundal et al. (2020). pVACtools: A computational toolkit to identify and visualize cancer neoantigens. Cancer Immunology Research.',
];
