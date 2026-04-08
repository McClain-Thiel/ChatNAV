import type { JobSummary } from '../api/types';
import JobCard from '../components/JobCard';

interface Props {
  jobs: JobSummary[];
}

export default function JobsPage({ jobs }: Props) {
  return (
    <div style={{
      maxWidth: 760,
      margin: '0 auto',
      padding: '40px 40px',
      height: '100%',
      overflowY: 'auto',
    }}>
      <div style={{
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'baseline',
        marginBottom: 24,
      }}>
        <h1 style={{
          fontFamily: "'EB Garamond', serif",
          fontSize: 28,
          fontWeight: 500,
          color: 'var(--text-primary)',
        }}>
          Jobs
        </h1>
        <button style={{
          fontFamily: "'EB Garamond', serif",
          fontSize: 14,
          color: 'var(--text-primary)',
          background: 'transparent',
          border: '1px solid var(--border-mid)',
          padding: '6px 14px',
          borderRadius: 2,
          cursor: 'pointer',
        }}>
          + New Job
        </button>
      </div>

      <div style={{ borderTop: '2px solid var(--border-heavy)', marginBottom: 2 }} />
      <div style={{ borderTop: '1px solid var(--border-mid)', marginBottom: 16 }} />

      {jobs.length === 0 ? (
        <div style={{
          textAlign: 'center',
          padding: '60px 0',
          fontFamily: "'Spectral', serif",
          fontStyle: 'italic',
          fontSize: 16,
          color: 'var(--text-meta)',
        }}>
          No jobs submitted yet. Click "+ New Job" to begin.
        </div>
      ) : (
        <div style={{
          border: '1px solid var(--border-light)',
          borderRadius: 3,
          overflow: 'hidden',
        }}>
          {jobs.map(job => (
            <JobCard key={job.id} job={job} />
          ))}
        </div>
      )}
    </div>
  );
}
