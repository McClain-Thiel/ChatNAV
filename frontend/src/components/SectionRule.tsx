interface Props {
  label?: string;
  style?: React.CSSProperties;
}

export default function SectionRule({ label, style }: Props) {
  return (
    <div style={{ marginBottom: 16, ...style }}>
      <div style={{ borderTop: '2px solid var(--border-heavy)' }} />
      {label && (
        <div style={{
          display: 'flex',
          justifyContent: 'center',
          marginTop: -8,
        }}>
          <span style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 8.5,
            letterSpacing: 2.5,
            textTransform: 'uppercase',
            color: 'var(--text-meta)',
            background: 'var(--bg-card)',
            padding: '0 12px',
          }}>
            {label}
          </span>
        </div>
      )}
      <div style={{ borderTop: '1px solid var(--border-mid)', marginTop: label ? 4 : 4 }} />
    </div>
  );
}
