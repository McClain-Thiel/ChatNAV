import { useState, useRef } from 'react';

interface Props {
  text: string;
  children: React.ReactNode;
  delay?: number;
}

export default function Tooltip({ text, children, delay = 600 }: Props) {
  const [visible, setVisible] = useState(false);
  const timerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const show = () => {
    timerRef.current = setTimeout(() => setVisible(true), delay);
  };

  const hide = () => {
    if (timerRef.current) clearTimeout(timerRef.current);
    setVisible(false);
  };

  return (
    <span
      onMouseEnter={show}
      onMouseLeave={hide}
      style={{ position: 'relative', cursor: 'help' }}
    >
      {children}
      {visible && (
        <span style={{
          position: 'absolute',
          bottom: '100%',
          left: '50%',
          transform: 'translateX(-50%)',
          marginBottom: 6,
          padding: '8px 12px',
          background: 'var(--text-primary)',
          color: '#faf7f2',
          fontFamily: "'Spectral', serif",
          fontSize: 12,
          lineHeight: 1.5,
          borderRadius: 3,
          whiteSpace: 'normal',
          width: 260,
          textAlign: 'left',
          fontWeight: 400,
          fontStyle: 'normal',
          letterSpacing: 0,
          textTransform: 'none',
          zIndex: 50,
          boxShadow: '0 4px 12px rgba(28,20,16,0.2)',
          animation: 'reveal 0.12s ease',
        }}>
          {text}
        </span>
      )}
    </span>
  );
}
