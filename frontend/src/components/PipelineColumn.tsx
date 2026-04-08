import type { JobDetail } from '../api/types';
import { getModuleSegments } from '../api/helpers';
import ProgressBar from './ProgressBar';
import ModuleRow from './ModuleRow';

interface Props {
  job: JobDetail;
  activeTab: 'pipeline' | 'results';
  onTabChange: (tab: 'pipeline' | 'results') => void;
  selectedModule: number;
  onModuleSelect: (n: number) => void;
}

export default function PipelineColumn({ job, activeTab, onTabChange, selectedModule, onModuleSelect }: Props) {
  const tumorType = (job.params.tumor_type as string) || '';

  return (
    <div style={{
      width: 258,
      flexShrink: 0,
      background: 'var(--bg-page)',
      borderRight: '1px solid var(--border-light)',
      display: 'flex',
      flexDirection: 'column',
      height: '100%',
      overflow: 'hidden',
    }}>
      {/* Header */}
      <div style={{ padding: '16px 14px 12px', borderBottom: '1px solid var(--border-light)' }}>
        <div style={{
          fontFamily: "'EB Garamond', serif",
          fontSize: 15,
          fontWeight: 600,
          color: 'var(--text-primary)',
          marginBottom: 4,
        }}>
          {job.patient_id}
        </div>
        <div style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 9,
          color: 'var(--text-muted)',
          marginBottom: 8,
        }}>
          {tumorType}
          {tumorType && ' · '}
          {new Date(job.created_at).toLocaleDateString('en-US', {
            day: 'numeric', month: 'long', year: 'numeric',
          })}
        </div>
        <ProgressBar segments={getModuleSegments(job.status, job.modules)} height={3} />

        {/* Tabs */}
        <div style={{ display: 'flex', gap: 16, marginTop: 12 }}>
          {(['results', 'pipeline'] as const).map(tab => (
            <button
              key={tab}
              onClick={() => onTabChange(tab)}
              style={{
                fontFamily: "'Inconsolata', monospace",
                fontSize: 9,
                letterSpacing: 2,
                textTransform: 'uppercase',
                background: 'none',
                border: 'none',
                padding: '4px 0',
                color: activeTab === tab ? 'var(--text-primary)' : 'var(--text-muted)',
                borderBottom: activeTab === tab ? '1.5px solid var(--text-primary)' : '1.5px solid transparent',
                cursor: 'pointer',
              }}
            >
              {tab}
            </button>
          ))}
        </div>
      </div>

      {/* Module list */}
      <div style={{ flex: 1, overflowY: 'auto' }}>
        {job.modules.map(mod => (
          <ModuleRow
            key={mod.number}
            module={mod}
            selected={selectedModule === mod.number}
            onClick={() => onModuleSelect(mod.number)}
          />
        ))}
      </div>
    </div>
  );
}
