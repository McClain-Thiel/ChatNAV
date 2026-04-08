import { BrowserRouter, Routes, Route } from 'react-router-dom';
import { useJobs } from './hooks/useJobs';
import Masthead from './components/Masthead';
import JobsPage from './pages/JobsPage';
import HomePage from './pages/HomePage';
import JobDetailPage from './pages/JobDetailPage';

export default function App() {
  const { jobs } = useJobs();

  return (
    <BrowserRouter>
      <div style={{ height: '100vh', display: 'flex', flexDirection: 'column' }}>
        <Masthead jobs={jobs} />
        <div style={{ marginTop: 56, flex: 1, overflow: 'hidden' }}>
          <Routes>
            <Route path="/" element={<HomePage />} />
            <Route path="/jobs" element={<JobsPage jobs={jobs} />} />
            <Route path="/jobs/:jobId" element={<JobDetailPage jobs={jobs} />} />
            <Route path="/jobs/:jobId/results" element={<JobDetailPage jobs={jobs} />} />
          </Routes>
        </div>
      </div>
    </BrowserRouter>
  );
}
