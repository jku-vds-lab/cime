export declare class CIMEBackend {
    readonly baseUrl: string;
    readonly fetchParams: Parameters<typeof fetch>[1];
    protected smiles_cache: {};
    protected smiles_highlight_cache: {};
    protected cache: {};
    constructor(baseUrl: string, fetchParams?: Parameters<typeof fetch>[1]);
    protected handleSmilesCache: (smiles: string, highlight?: boolean) => any;
    protected setSmilesCache: (smiles: any, highlight: boolean, data: any) => void;
    protected async_cache: (cached_data: any) => Promise<any>;
    handleCache: (key: any) => any;
    setCache: (key: any, value: any) => void;
    handleErrors: (response: any) => any;
    handleJSONErrors: (data: any) => any;
    deleteFile: (filename: any) => Promise<{
        deleted: boolean;
    }>;
    getFiles: () => Promise<{
        name: string;
        id: number;
    }[]>;
    getDifferenceHighlight: (smilesA: any, smilesB: any, controller: any) => Promise<any>;
    getStructureFromSmiles: (id: string | number, smiles: string, highlight: boolean, controller: any) => Promise<any>;
    getStructuresFromSmilesList: (formData: FormData, controller?: any) => Promise<any>;
    getMCSFromSmilesList: (formData: FormData, controller?: any) => Promise<any>;
    getSubstructureCount: (smiles_list: any, filter: any) => Promise<any>;
    upload_sdf_file: (file: any, controller?: any) => Promise<{
        name: string;
        id: number;
    }>;
    getRepresentationList: (refresh: boolean, id: string | number, controller: AbortController) => Promise<any>;
    calculateHDBScanClusters: (X: any, min_cluster_size: any, min_cluster_samples: any, allow_single_cluster: any) => Promise<any>;
    calculateScores: (id: string, cime_ids: string[], currentRep: string) => Promise<any>;
}
export declare const CIMEBackendFromEnv: CIMEBackend;
//# sourceMappingURL=CIMEBackend.d.ts.map