import { useNavigate } from 'react-router-dom';
import type { JobSummary } from '../api/types';
import StatusBadge from './StatusBadge';
import ProgressBar from './ProgressBar';
import { getModuleSegments } from '../api/helpers';

interface Props {
  job: JobSummary;
  compact?: boolean;
  active?: boolean;
}

export default function JobCard({ job, compact = false, active = false }: Props) {
  const navigate = useNavigate();
  const tumorType = (job.params.tumor_type as string) || '';

  return (
    <div
      onClick={() => navigate(`/jobs/${job.id}`)}
      style={{
        padding: compact ? '10px 14px' : '14px 16px',
        borderLeft: active ? '2px solid var(--border-heavy)' : '2px solid transparent',
        borderBottom: '1px solid var(--border-faint)',
        background: active ? 'var(--bg-active)' : 'transparent',
        cursor: 'pointer',
        transition: 'background 0.1s',
      }}
      onMouseEnter={e => {
        if (!active) e.currentTarget.style.background = 'var(--bg-hover)';
      }}
      onMouseLeave={e => {
        if (!active) e.currentTarget.style.background = 'transparent';
      }}
    >
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: 4 }}>
        <span style={{
          fontFamily: "'EB Garamond', serif",
          fontSize: compact ? 14 : 16,
          fontWeight: 600,
          color: 'var(--text-primary)',
        }}>
          {job.patient_id}
        </span>
        <StatusBadge status={job.status} />
      </div>
      {tumorType && (
        <div style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: compact ? 9 : 10,
          color: 'var(--text-meta)',
          marginBottom: 4,
        }}>
          {tumorType}
        </div>
      )}
      <div style={{
        fontFamily: "'Inconsolata', monospace",
        fontSize: 9,
        color: 'var(--text-faint)',
        marginBottom: 6,
      }}>
        {new Date(job.created_at).toLocaleDateString('en-US', {
          day: 'numeric', month: 'long', year: 'numeric',
        })}
      </div>
      <ProgressBar segments={getModuleSegments(job.status)} height={compact ? 2 : 3} />
    </div>
  );
}
