interface Props {
  /** Array of segment statuses: 'complete', 'running', 'queued' */
  segments: ('complete' | 'running' | 'queued')[];
  height?: number;
}

const SEGMENT_COLORS = {
  complete: 'var(--status-complete)',
  running: 'var(--status-running)',
  queued: 'var(--border-faint)',
};

export default function ProgressBar({ segments, height = 3 }: Props) {
  return (
    <div style={{ display: 'flex', gap: 2 }}>
      {segments.map((seg, i) => (
        <div
          key={i}
          style={{
            flex: 1,
            height,
            borderRadius: 1,
            background: SEGMENT_COLORS[seg],
          }}
        />
      ))}
    </div>
  );
}
