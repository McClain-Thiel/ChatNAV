import type { JobSummary, JobDetail, OutputFile } from './types';

const API_BASE = '/api/v1';

async function fetchJSON<T>(path: string): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`);
  if (!res.ok) {
    throw new Error(`API error: ${res.status} ${res.statusText}`);
  }
  return res.json();
}

export async function getJobs(): Promise<JobSummary[]> {
  const data = await fetchJSON<{ jobs: JobSummary[]; total: number }>('/jobs?limit=100');
  return data.jobs;
}

export async function getJob(jobId: string): Promise<JobDetail> {
  return fetchJSON<JobDetail>(`/jobs/${jobId}`);
}

export async function getJobOutputs(jobId: string): Promise<OutputFile[]> {
  const data = await fetchJSON<{ files: OutputFile[] }>(`/jobs/${jobId}/outputs`);
  return data.files;
}

export async function getJobLog(jobId: string, tail = 0): Promise<string> {
  const data = await fetchJSON<{ log: string }>(`/jobs/${jobId}/log?tail=${tail}`);
  return data.log;
}

export function getDownloadUrl(jobId: string, filePath: string): string {
  return `${API_BASE}/jobs/${jobId}/outputs/${filePath}`;
}

export async function cancelJob(jobId: string): Promise<void> {
  const res = await fetch(`${API_BASE}/jobs/${jobId}/cancel`, { method: 'POST' });
  if (!res.ok) {
    throw new Error(`Cancel failed: ${res.status}`);
  }
}

export interface FilePreview {
  type: 'table' | 'json' | 'fasta' | 'text';
  columns?: string[];
  rows?: string[][];
  total_rows?: number;
  data?: Record<string, unknown>;
  sequences?: { header: string; sequence: string }[];
  content?: string;
}

export async function getFilePreview(jobId: string, filePath: string): Promise<FilePreview> {
  return fetchJSON<FilePreview>(`/jobs/${jobId}/preview/${filePath}`);
}

export async function submitJob(formData: FormData): Promise<JobSummary> {
  const res = await fetch(`${API_BASE}/jobs`, {
    method: 'POST',
    body: formData,
  });
  if (!res.ok) {
    throw new Error(`Submit failed: ${res.status}`);
  }
  return res.json();
}
