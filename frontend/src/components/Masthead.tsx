import { NavLink, Link } from 'react-router-dom';
import type { JobSummary } from '../api/types';

interface Props {
  jobs: JobSummary[];
}

export default function Masthead({ jobs }: Props) {
  const hasRunning = jobs.some(j => j.status === 'running');

  return (
    <header style={styles.header}>
      <div style={styles.topRule} />
      <div style={styles.inner}>
        <Link to="/" style={{ ...styles.brand, textDecoration: 'none' }}>
          <span style={styles.wordmark}>ChatNAV</span>
          <span style={styles.subtitle}>PERSONALISED CANCER IMMUNOTHERAPY</span>
        </Link>
        <nav style={styles.nav}>
          {[
            { to: '/jobs', label: 'Jobs' },
            { to: '/reference', label: 'Reference' },
            { to: '/settings', label: 'Settings' },
          ].map(link => (
            <NavLink
              key={link.to}
              to={link.to}
              style={({ isActive }) => ({
                ...styles.navLink,
                fontWeight: isActive ? 600 : 400,
                color: isActive ? 'var(--text-primary)' : 'var(--text-muted)',
                borderBottom: isActive ? '1.5px solid var(--border-heavy)' : '1.5px solid transparent',
              })}
            >
              {link.label}
            </NavLink>
          ))}
        </nav>
        <div style={styles.right}>
          {hasRunning && (
            <span style={styles.statusPill}>
              <span className="animate-blink" style={styles.statusDot}>●</span>
              {' '}Running
            </span>
          )}
          <div style={styles.avatar}>MT</div>
        </div>
      </div>
      <div style={styles.bottomRule} />
    </header>
  );
}

const styles: Record<string, React.CSSProperties> = {
  header: {
    position: 'fixed',
    top: 0,
    left: 0,
    right: 0,
    height: 56,
    zIndex: 100,
    background: 'var(--bg-page)',
  },
  topRule: {
    height: 3,
    background: 'var(--border-heavy)',
  },
  inner: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    height: 52,
    padding: '0 24px',
  },
  brand: {
    display: 'flex',
    alignItems: 'baseline',
    gap: 12,
  },
  wordmark: {
    fontFamily: "'EB Garamond', serif",
    fontSize: 20,
    fontWeight: 600,
    color: 'var(--text-primary)',
  },
  subtitle: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 8.5,
    letterSpacing: 2.5,
    textTransform: 'uppercase' as const,
    color: 'var(--text-muted)',
  },
  nav: {
    display: 'flex',
    gap: 28,
  },
  navLink: {
    fontFamily: "'EB Garamond', serif",
    fontSize: 15,
    paddingBottom: 2,
    textDecoration: 'none',
  },
  right: {
    display: 'flex',
    alignItems: 'center',
    gap: 14,
  },
  statusPill: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 9,
    color: 'var(--status-running)',
  },
  statusDot: {
    color: 'var(--status-running)',
  },
  avatar: {
    width: 26,
    height: 26,
    borderRadius: '50%',
    background: 'linear-gradient(135deg, var(--border-mid), var(--text-muted))',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    fontFamily: "'Inconsolata', monospace",
    fontSize: 9,
    fontWeight: 500,
    color: '#fff',
  },
  bottomRule: {
    height: 1,
    background: 'var(--border-light)',
  },
};
