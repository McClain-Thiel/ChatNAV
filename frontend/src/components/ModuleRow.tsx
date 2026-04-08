import type { ModuleState } from '../api/types';
import { formatDuration, getEpitopeColor } from '../api/helpers';

interface Props {
  module: ModuleState;
  selected: boolean;
  onClick: () => void;
}

export default function ModuleRow({ module, selected, onClick }: Props) {
  const isComplete = module.status === 'complete';
  const isRunning = module.status === 'running';
  const isQueued = module.status === 'queued';

  return (
    <div
      onClick={onClick}
      style={{
        display: 'flex',
        alignItems: 'flex-start',
        gap: 10,
        padding: '10px 14px',
        borderLeft: selected ? `2px solid ${getEpitopeColor(module.number)}` : '2px solid transparent',
        borderBottom: '1px solid var(--border-faint)',
        background: selected ? 'var(--bg-active)' : 'transparent',
        cursor: 'pointer',
        transition: 'background 0.1s',
      }}
      onMouseEnter={e => {
        if (!selected) e.currentTarget.style.background = 'var(--bg-hover)';
      }}
      onMouseLeave={e => {
        if (!selected) e.currentTarget.style.background = 'transparent';
      }}
    >
      {/* Status dot */}
      <div style={{
        width: 10,
        height: 10,
        borderRadius: '50%',
        marginTop: 4,
        flexShrink: 0,
        ...(isComplete && { background: 'var(--status-complete)' }),
        ...(isRunning && { background: 'var(--status-running)' }),
        ...(isQueued && { border: '1.5px solid var(--border-light)', background: 'transparent' }),
        ...(module.status === 'failed' && { background: 'var(--status-error)' }),
      }}>
        {isRunning && <div className="animate-blink" style={{ width: '100%', height: '100%', borderRadius: '50%', background: 'var(--status-running)' }} />}
      </div>

      {/* Step number */}
      <span style={{
        fontFamily: "'Inconsolata', monospace",
        fontSize: 9,
        color: 'var(--text-faint)',
        minWidth: 16,
        flexShrink: 0,
        marginTop: 3,
      }}>
        {String(module.number).padStart(2, '0')}
      </span>

      {/* Title + subtitle */}
      <div style={{ flex: 1, minWidth: 0 }}>
        <div style={{
          fontFamily: "'EB Garamond', serif",
          fontSize: 13.5,
          fontWeight: isQueued ? 400 : 600,
          color: isQueued ? 'var(--text-faint)' : 'var(--text-primary)',
          lineHeight: 1.3,
        }}>
          {module.name}
        </div>
        <div style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 9,
          color: 'var(--text-faint)',
          marginTop: 2,
        }}>
          {module.subtitle}
        </div>
        {isRunning && module.progress_pct != null && (
          <div style={{
            marginTop: 6,
            height: 2,
            background: 'var(--border-light)',
            borderRadius: 1,
          }}>
            <div style={{
              height: '100%',
              width: `${module.progress_pct}%`,
              background: 'var(--status-running)',
              borderRadius: 1,
              transition: 'width 0.3s ease',
            }} />
          </div>
        )}
      </div>

      {/* Duration */}
      {isComplete && module.duration_seconds != null && (
        <span style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 8.5,
          color: 'var(--status-complete)',
          flexShrink: 0,
          marginTop: 3,
        }}>
          {formatDuration(module.duration_seconds)}
        </span>
      )}
    </div>
  );
}
