import type { JobSummary } from '../api/types';
import JobCard from './JobCard';

interface Props {
  jobs: JobSummary[];
  activeJobId: string | undefined;
}

export default function JobSidebar({ jobs, activeJobId }: Props) {
  return (
    <div style={{
      width: 220,
      flexShrink: 0,
      background: 'var(--bg-sidebar)',
      borderRight: '1px solid var(--border-light)',
      display: 'flex',
      flexDirection: 'column',
      height: '100%',
      overflow: 'hidden',
    }}>
      {/* Header */}
      <div style={{
        padding: '14px 14px',
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        letterSpacing: 2.5,
        textTransform: 'uppercase',
        color: 'var(--text-muted)',
        borderBottom: '1px solid var(--border-light)',
      }}>
        Jobs
      </div>

      {/* Job list */}
      <div style={{ flex: 1, overflowY: 'auto' }}>
        {jobs.map(job => (
          <JobCard
            key={job.id}
            job={job}
            compact
            active={job.id === activeJobId}
          />
        ))}
      </div>

      {/* New job button */}
      <button style={{
        margin: 10,
        padding: '10px 0',
        fontFamily: "'Inconsolata', monospace",
        fontSize: 9.5,
        color: 'var(--text-muted)',
        background: 'transparent',
        border: '1px dashed var(--border-light)',
        borderRadius: 2,
        cursor: 'pointer',
      }}>
        + New Job
      </button>
    </div>
  );
}
