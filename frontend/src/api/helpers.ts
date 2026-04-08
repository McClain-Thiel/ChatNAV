import type { JobStatus, ModuleState } from './types';

const EPITOPE_COLORS = [
  'var(--ep-01)', 'var(--ep-02)', 'var(--ep-03)', 'var(--ep-04)',
  'var(--ep-05)', 'var(--ep-06)', 'var(--ep-07)', 'var(--ep-08)',
  'var(--ep-09)', 'var(--ep-10)', 'var(--ep-11)', 'var(--ep-12)',
  'var(--ep-13)', 'var(--ep-14)', 'var(--ep-15)', 'var(--ep-16)',
  'var(--ep-17)', 'var(--ep-18)', 'var(--ep-19)', 'var(--ep-20)',
];

export function getEpitopeColor(rank: number): string {
  return EPITOPE_COLORS[(rank - 1) % EPITOPE_COLORS.length];
}

export function getScoreColor(score: number): string {
  if (score >= 0.70) return 'var(--score-high)';
  if (score >= 0.50) return 'var(--score-mid)';
  return 'var(--score-low)';
}

export function formatDuration(seconds: number): string {
  if (seconds < 60) return `${Math.round(seconds)}s`;
  if (seconds < 3600) return `${Math.round(seconds / 60)}m`;
  const h = Math.floor(seconds / 3600);
  const m = Math.round((seconds % 3600) / 60);
  return m > 0 ? `${h}h ${m}m` : `${h}h`;
}

/** Generate 11 progress bar segments from job status or module states */
export function getModuleSegments(
  status: JobStatus,
  modules?: ModuleState[],
): ('complete' | 'running' | 'queued')[] {
  if (modules && modules.length > 0) {
    return modules.map(m => {
      if (m.status === 'complete') return 'complete';
      if (m.status === 'running') return 'running';
      return 'queued';
    });
  }

  // Fallback based on job status
  if (status === 'completed') return Array(11).fill('complete');
  if (status === 'failed') return Array(11).fill('queued');
  return Array(11).fill('queued');
}

/** Compute P(≥1 neoantigen active) */
export function computeProbability(ppvPerEpitope: number, nEpitopes: number): number {
  return 1 - Math.pow(1 - ppvPerEpitope, nEpitopes);
}

/** Get TCR-facing positions for a peptide of given length */
export function getTCRFacingPositions(length: number): number[] {
  // Positions 3 through length-3 (0-indexed) are TCR-facing
  const positions: number[] = [];
  for (let i = 3; i <= length - 3; i++) {
    positions.push(i);
  }
  return positions;
}

/** Get the mutated residue position (center of peptide) */
export function getMutatedPosition(length: number): number {
  return Math.floor(length / 2);
}
