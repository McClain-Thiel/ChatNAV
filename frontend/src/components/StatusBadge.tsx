import type { JobStatus, ModuleStatus } from '../api/types';

type Status = JobStatus | ModuleStatus;

const STATUS_CONFIG: Record<string, { color: string; bg: string; label: string }> = {
  completed: { color: 'var(--status-complete)', bg: 'var(--status-complete-bg)', label: 'COMPLETE' },
  complete: { color: 'var(--status-complete)', bg: 'var(--status-complete-bg)', label: 'COMPLETE' },
  running: { color: 'var(--status-running)', bg: 'var(--status-running-bg)', label: 'RUNNING' },
  pending: { color: 'var(--text-faint)', bg: 'transparent', label: 'QUEUED' },
  queued: { color: 'var(--text-faint)', bg: 'transparent', label: 'QUEUED' },
  failed: { color: 'var(--status-error)', bg: 'var(--status-error-bg)', label: 'FAILED' },
  cancelled: { color: 'var(--text-faint)', bg: 'transparent', label: 'CANCELLED' },
};

interface Props {
  status: Status;
}

export default function StatusBadge({ status }: Props) {
  const config = STATUS_CONFIG[status] ?? STATUS_CONFIG.pending;

  return (
    <span
      style={{
        fontFamily: "'Inconsolata', monospace",
        fontSize: 9,
        letterSpacing: 1.5,
        textTransform: 'uppercase',
        color: config.color,
        background: config.bg,
        padding: '2px 6px',
        borderRadius: 2,
      }}
    >
      {status === 'running' && (
        <span className="animate-blink" style={{ marginRight: 4 }}>●</span>
      )}
      {status === 'completed' || status === 'complete' ? '✓ ' : ''}
      {config.label}
    </span>
  );
}
