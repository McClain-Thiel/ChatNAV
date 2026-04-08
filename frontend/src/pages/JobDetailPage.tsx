import { useState } from 'react';
import { useParams } from 'react-router-dom';
import type { JobSummary } from '../api/types';
import { useJob } from '../hooks/useJob';
import { HoverContext } from '../context/HoverContext';
import JobSidebar from '../components/JobSidebar';
import PipelineColumn from '../components/PipelineColumn';
import PipelineDetail from '../components/PipelineDetail';
import ResultsDetail from '../components/ResultsDetail';

interface Props {
  jobs: JobSummary[];
}

export default function JobDetailPage({ jobs }: Props) {
  const { jobId } = useParams<{ jobId: string }>();
  const { job, loading, error } = useJob(jobId);
  const [activeTab, setActiveTab] = useState<'pipeline' | 'results'>('results');
  const [selectedModule, setSelectedModule] = useState(1);
  const [hoveredRank, setHoveredRank] = useState<number | null>(null);

  if (loading && !job) {
    return (
      <div style={{ display: 'flex', height: '100%' }}>
        <JobSidebar jobs={jobs} activeJobId={jobId} />
        <div style={{
          flex: 1,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontFamily: "'Spectral', serif",
          fontStyle: 'italic',
          color: 'var(--text-meta)',
        }}>
          Loading...
        </div>
      </div>
    );
  }

  if (error || !job) {
    return (
      <div style={{ display: 'flex', height: '100%' }}>
        <JobSidebar jobs={jobs} activeJobId={jobId} />
        <div style={{
          flex: 1,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontFamily: "'Spectral', serif",
          color: 'var(--status-error)',
        }}>
          {error || 'Job not found'}
        </div>
      </div>
    );
  }

  const currentModule = job.modules.find(m => m.number === selectedModule) ?? job.modules[0];

  return (
    <HoverContext.Provider value={{ hoveredRank, setHoveredRank }}>
      <div style={{ display: 'flex', height: '100%' }}>
        <JobSidebar jobs={jobs} activeJobId={jobId} />
        <PipelineColumn
          job={job}
          activeTab={activeTab}
          onTabChange={setActiveTab}
          selectedModule={selectedModule}
          onModuleSelect={setSelectedModule}
        />
        <div style={{
          flex: 1,
          overflowY: 'auto',
          background: 'var(--bg-card)',
        }}>
          {activeTab === 'pipeline' && currentModule && (
            <PipelineDetail
              module={currentModule}
              totalModules={job.modules.length}
              jobId={job.id}
            />
          )}
          {activeTab === 'results' && (
            <ResultsDetail job={job} />
          )}
        </div>
      </div>
    </HoverContext.Provider>
  );
}
